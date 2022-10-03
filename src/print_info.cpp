#include "print_info.h"

void summarize_sing_sra(SRA sra) {
  std::cout << sra.get_accession() << " - "
            << sra.get_org_name()  << " " 
            << sra.get_tax_id()    << std::endl;
}

void summarize_all_sras(std::vector<SRA> sras) {
  for (auto sra : sras) {
    summarize_sing_sra(sra);
  }
}

