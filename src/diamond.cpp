#include "diamond.h"

void makeDb(std::string pathProtRef, std::string outDir) {
  fs::path pathProtRefFile(pathProtRef.c_str());
  std::string protRefStr(pathProtRefFile.stem().c_str());
  fs::path outDbFile((outDir + "/" + protRefStr).c_str());
  int result;
  if (fs::exists(outDbFile)) {
    std::cout << "BLAST database found for: " + protRefStr << std::endl;
    return;
  }
  std::string makeDbCmd = DMAKEDB + " --in " + pathProtRef + " --db " +
                          outDbFile.c_str();
  result = system(makeDbCmd.c_str());
  if (WIFSIGNALED(result)) {
    std::cout << "Existed with signal " << WTERMSIG(result) << std::endl;
    exit(1);
  } 
}

void blastxDiam(transcript transcripts, std::string blastDb,
                std::string threads, std::string outDir) {
  fs::path pathTrans = transcripts.get_trans_path_trinity();
  std::string pathTransStr(pathTrans.c_str());
  std::string transStr(pathTrans.stem().c_str());

  fs::path outBlastxFile((outDir + "/" + transStr + ".blastx").c_str());
  std::string outBlastxStr(outBlastxFile.c_str());
  int result;
  if (fs::exists(outBlastxFile)) {
    std::cout << "Diamond BLASTX output found for: " << transStr << std::endl;
    return;
  }
  std::string blastxDiamCmd = DBLASTX + " --db " + blastDb + " --query " + pathTransStr +
                              " --evalue " + "0.01" +
                              " --outfmt " + "6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend start send evalue bitscore" +
                              " --out " + outBlastxStr +
                              " --threads " + threads +
                              " --max_target_seqs 100";
  result = system(blastxDiamCmd.c_str());
  if (WIFSIGNALED(result)) {
    std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
    exit(1);
  }
}

void blastpDiam(std::string pepFilePath, std::string blastDb,
                std::string threads, std::string outFile) {
  fs::path outBlastpFile(outFile.c_str());
  int result;
  if (fs::exists(outBlastpFile)) {
    std::cout << "Diamond BLASTP output found for: " << pepFilePath << std::endl;
    return;
  }
  std::string  blastpDiamCmd = DBLASTP + " --query " + pepFilePath + " --db " + blastDb +
                               " --max_target_seqs 1 --outfmt 6 --evalue 10 " +
                               " --threads " + threads + " --out " + outFile;
  result = system(blastpDiamCmd.c_str());
  if (WIFSIGNALED(result)) {
    std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
    exit(1);
  }
}
