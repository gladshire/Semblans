#include "sra_toolkit.h"

namespace fs = boost::filesystem;


void prefetch_sra(std::vector<SRA> sras) {
  std::string outDir(sras[0].get_sra_path_raw().first.parent_path().native().c_str());
  std::string prefetchFlag = " --max-size u -p -o ";
  for (auto sra : sras) {
    system((PATH_PREFETCH + " " + sra.get_accession() + prefetchFlag + outDir).c_str());
  }
}

void fasterq_sra(std::vector<SRA> sras, std::string threads) {
  std::string prefetchDir(sras[0].get_sra_path_raw().first.parent_path().native().c_str());
  std::string outFile;
  std::string fasterqFlag = " -p -e " + threads + " " + prefetchDir;
  for (auto sra : sras) {
    outFile = sra.make_file_str();
    system((PATH_FASTERQ + " " + fasterqFlag + sra.get_accession() + " -o " + prefetchDir + outFile).c_str());
  }
}

std::vector<SRA> get_sras(const INI_MAP &iniFile) {
  std::vector<SRA> sras;
  for (auto sra : iniFile.at("SRA accessions")) {
    sras.push_back(SRA(sra.first));
  }
  return sras;
}


int main() {
  extern INI_MAP cfgIni;
  std::vector<SRA> sras = get_sras(cfgIni);
  fs::path outdir = cfgIni["General"]["output_directory"];
  prefetch_sra(sras);
  fasterq_sra(sras, "2");
  return 0;
}
