#include "print_info.h"

void summarize_sing_sra(SRA sra) {
  std::string filePrefix1 = sra.get_file_prefix().first;
  std::string filePrefix2 = sra.get_file_prefix().second;
  if (sra.is_paired()) {
    std::cout << "  Paired-end run:" << std::endl;
    std::cout << "  " << filePrefix1 << std::endl;
    std::cout << "  " << filePrefix2 << "\n" << std::endl;
  }
  else {
    std::cout << "  Single-end run:" << std::endl;
    std::cout << "  " << filePrefix1 << "\n" << std::endl;
  }
}

void summarize_all_sras(const std::vector<SRA> & sras) {
  for (auto sra : sras) {
    summarize_sing_sra(sra);
  }
}

