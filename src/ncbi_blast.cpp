#include "ncbi_blast.h"

void makeBlastDb(std::string pathProtRef, std::string outDir,
                 bool dispOutput, std::string logFile) {
  fs::path pathProtRefFile(pathProtRef.c_str());
  std::string protRefStr(pathProtRefFile.stem().c_str()); 
  fs::path outDbFile((outDir + "/" + protRefStr).c_str());
  std::string printOut;
  if (dispOutput) {
    printOut = " 2&>1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  if (fs::exists(outDbFile)) {
    logOutput("BLAST database found for: " + protRefStr, logFile);
    return;
  }
  std::string makeblastdbCmd = MAKEBLASTDB + " -in " + pathProtRef + " -dbtype prot -out " +
                               outDbFile.c_str() + printOut;
  result = system(makeblastdbCmd.c_str());
  checkExitSignal(result, logFile);
}


void blastx(std::string transIn, std::string blastDb,
            std::string maxEvalue, std::string maxTargetSeqs,
            std::string threads, std::string outDir,
            bool dispOutput, std::string logFile) {
  fs::path pathTrans(transIn.c_str());
  std::string pathTransStr(pathTrans.c_str());
  std::string transStr(pathTrans.stem().c_str());

  fs::path outBlastxFile((outDir + "/" + transStr + ".blastx").c_str());
  std::string outBlastxStr(outBlastxFile.c_str());
  std::string printOut;
  if (dispOutput) {
    printOut = " 2&>1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  if (fs::exists(outBlastxFile)) {
    logOutput("BLASTX output found for: " + transStr, logFile);
    return;
  }
  std::string blastxCmd = BLASTX + " -db " + blastDb + " -query " + pathTransStr +
                          " -evalue " + maxEvalue + " -max_target_seqs " + maxTargetSeqs +
                          " -outfmt " + "\"6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore\"" +
                          " -out " + outBlastxStr +
                          " -num_threads " + threads +
                          " -max_target_seqs 100" + printOut;
  result = system(blastxCmd.c_str());
  checkExitSignal(result, logFile);
}


void blastp(std::string pepFilePath, std::string blastDb,
            std::string maxEvalue, std::string maxTargetSeqs,
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
    logOutput("BLASTP output found for: " + pepFilePath, logFile);
    return;
  }
  std::string blastpCmd = BLASTP + " -query " + pepFilePath + " -db " + blastDb +
                          " -evalue " + maxEvalue + " -max_target_seqs " + maxTargetSeqs + 
                          " -outfmt 6 " + "-num_threads " + threads + " > " + outFile + printOut;
  result = system(blastpCmd.c_str());
  checkExitSignal(result, logFile);
}
