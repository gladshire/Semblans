#include "trimm_wrap.h"

void run_trimmomatic(const std::vector<SRA> & sras, std::string threads) {
  std::cout << "\nTrimming adapter sequences for:\n" << std::endl;
  summarize_all_sras(sras);

  std::string inFile1;
  std::string inFile2;
  std::string outFileP1;
  std::string outFileP2;
  std::string outFileU1;
  std::string outFileU2;
  std::string outFile;
  int result;
  for (auto sra : sras) {
    std::string outDir(sra.get_sra_path_trim_u().first.parent_path().c_str());
    std::string trimmCmd("java -jar " + PATH_TRIMM);
    std::string trimmFlags("-threads " + threads + " " + "ILLUMINACLIP:" + TRUSEQ_ALL +
                           ":2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25");
    inFile1 = sra.get_sra_path_corr_fix().first.c_str();
    if (sra.is_paired()) {
      inFile2 = sra.get_sra_path_corr_fix().second.c_str();
      outFileP1 = sra.get_sra_path_trim_p().first.c_str();
      outFileP2 = sra.get_sra_path_trim_p().second.c_str();
      outFileU1 = sra.get_sra_path_trim_u().first.c_str();
      outFileU2 = sra.get_sra_path_trim_u().second.c_str();
      if (fs::exists(sra.get_sra_path_trim_p().first) &&
          fs::exists(sra.get_sra_path_trim_p().second)) {
        std::cout << "Trimmed version found for: " << sra.get_accession() << std::endl;
        continue;
      }
      result = system((trimmCmd + " PE " + inFile1 + " " + inFile2 + " " + outFileP1 + " " + outFileU1 +
                        " " + outFileP2 + " " + outFileU2 + " " + trimmFlags).c_str());
    }
    else {
      outFile = sra.get_sra_path_trim_u().first.c_str();
      if (fs::exists(sra.get_sra_path_trim_p().first)) {
        std::cout << "Trimmed version found for " << sra.get_accession() << std::endl;
        continue;
      result = system((trimmCmd + " SE " + inFile1 + " " + outFile + " " + trimmFlags).c_str());
      }
    }
    if (WIFSIGNALED(result)) {
      std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
      exit(1);
    }
  }
}
