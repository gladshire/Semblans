#include "rem_overrep.h"

std::vector<std::string> get_overrep_seqs(SRA sra) {
  std::string inFile1Str(std::string(sra.get_fastqc_dir().first.c_str()) + "/" +
                         sra.make_file_str() + "_fastqc.html");
  if (sra.is_paired()) {
    std::string inFile2Str(std::string(sra.get_fastqc_dir().second.c_str()) + "/" +
                           sra.make_file_str() + "_fastqc.html");
  }
  // Return over-represented sequences in vector form for SRA file, single or paired
}

