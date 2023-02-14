#include "transdecoder_wrap.h"

int MIN_FASTA = 1000;

bool fasta_ok(std::string fastaFile, uintmax_t ram_b) {
  fs::path fastaFilePath(fastaFile.c_str());
  if (!fs::exists(fastaFilePath)) {
    return false;
  }
  uintmax_t numBytesTrans = fs::file_size(fastaFilePath);
  uintmax_t lenHashTable = numBytesTrans / 160;
  uintmax_t fasta_count = 0;
  seqHash fastaHash(lenHashTable, fastaFilePath, ram_b);
  fasta_count = fastaHash.getSize();
  if (fasta_count >= MIN_FASTA) {
    return true;
  }
  else {
    return false;
  }
}

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

void run_transdecoder(transcript trans, std::string threads, uintmax_t ram_b,
                      std::string dbPath, std::string outDir) {
  std::string transFilePath(trans.get_trans_path_largest().c_str());
  fs::path cdsFilePath = trans.get_trans_path_cds().c_str();
  fs::path pepFilePath = trans.get_trans_path_prot().c_str();

  std::string blastpout = trans.make_file_str() + ".blastp.outfmt6";
  if (fasta_ok(std::string(cdsFilePath.c_str()), ram_b) &&
      fasta_ok(std::string(pepFilePath.c_str()), ram_b)) {
    std::cout << "Skipping TransDecoder" << std::endl;
  }
  else {
    std::string transFileName = std::string(fs::path(transFilePath.c_str()).filename().c_str());
    fs::path allpep(std::string(outDir + transFileName +
                    ".transdecoder_dir/longest_orfs.pep").c_str());
    if (fs::exists(allpep)) {
      std::cout << "Skipping search for long orfs" << std::endl;
    }
    else {
      // Only operates on un-stranded. Implement stranded later
      std::string tdLongOrfs_cmd = PATH_TRANSD_LONGORFS + " -t " + transFilePath + " -O " +
                                   std::string(allpep.parent_path().c_str());
      system(tdLongOrfs_cmd.c_str());
    }
    if (blastpout_ok(outDir + blastpout)) {
      std::cout << "Skipping blastp" << std::endl;
    }
    else {
      // Update blastp path. Add wrapper
      std::string blastp_cmd = std::string((dl::program_location().parent_path() /
                                            fs::path(std::string("../lib/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/blastp"))).c_str()) +
                               " -query " + std::string(allpep.c_str()) +
                               " -db " + dbPath + " -max_target_seqs 1 -outfmt 6 -evalue 10 " +
                               "-num_threads " + threads + " > " + outDir + blastpout;
      system(blastp_cmd.c_str());
    }
    if (fasta_ok(std::string(cdsFilePath.c_str()), ram_b) &&
        fasta_ok(std::string(pepFilePath.c_str()), ram_b)) {
      std::cout << "Skip finding final CDS and PEP" << std::endl;
    }
    else {
      std::string tdPredict_cmd = PATH_TRANSD_PREDICT + " -t " + transFilePath +
                                  " --retain_blastp_hits " + outDir + blastpout + " --cpu " +
                                  threads;
      system(tdPredict_cmd.c_str());
    }
  }
}
