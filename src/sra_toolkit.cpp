#include "sra_toolkit.h"

namespace fs = boost::filesystem;


void prefetch_sra(std::vector<SRA> sras) {
  std::string outDir(sras[0].get_sra_path_raw().first.parent_path().native().c_str());
  std::string prefetchFlag = " --max-size u -p -O ";
  for (auto sra : sras) {
    system((PATH_PREFETCH + " " + sra.get_accession() + prefetchFlag + outDir).c_str());
  }
}

void fasterq_sra(std::vector<SRA> sras, std::string threads) {
  std::string prefetchDir(sras[0].get_sra_path_raw().first.parent_path().native().c_str());
  std::string outFile;
  std::string fasterqFlag = " -p -e " + threads + " -t " + prefetchDir + " " + prefetchDir;
  for (auto sra : sras) {
    outFile = sra.make_file_str();
    system((PATH_FASTERQ + " " + fasterqFlag + "/" + sra.get_accession() +
            " -o " + prefetchDir + "/" + outFile).c_str());
  }
}

