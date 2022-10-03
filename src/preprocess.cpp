#include "preprocess.h"

std::vector<std::string> stepDirs = {"00-Raw_reads/", "01-Quality_analysis_1/",
                                     "02-Error_correction/", "03-Trimming/",
                                     "04-Filtering/"};


std::vector<SRA> get_sras(const INI_MAP &iniFile) {
  std::vector<SRA> sras;
  for (auto sra : iniFile.at("SRA accessions")) {
    sras.push_back(SRA(sra.first, iniFile));
  }
  return sras;
}

void retrieve_sra_data(std::vector<SRA> sras, std::string threads) {
  std::cout << "Retrieving SRA runs for:\n" << std::endl;
  summarize_all_sras(sras);

  prefetch_sra(sras);
  fasterq_sra(sras, threads);
}

void make_proj_space(const INI_MAP &iniFile) {
  std::string projDir = iniFile.at("General").at("output_directory") +
                        iniFile.at("General").at("project_name") + "/";
  extern std::vector<std::string> stepDirs;
  system(("mkdir " + projDir).c_str());
  for (auto dir : stepDirs) {
    system(("mkdir " + projDir + dir).c_str());
  }
}

void print_help() {
  std::cout << "\n" << "NAME_OF_PROGRAM" << " - "
            << "A tool for bulk assemblies of de novo transcriptome data" << std::endl;
  std::cout << "\n" << "COMMAND STRUCTURE" << std::endl;
  std::cout << "\n" << "preprocess PATH/TO/CONFIG.INI num_threads" << std::endl;
}
int main(int argc, char * argv[]) {
  if (argc != 3) {
    print_help();
    return 0;
  }
  else {
    INI_MAP cfgIni = make_ini_map(argv[1]);
    std::string threads = argv[2];
    std::vector<SRA> sras = get_sras(cfgIni);
    
    make_proj_space(cfgIni);
    retrieve_sra_data(sras, threads);
    run_fastqc(sras, threads);
    run_rcorr(sras, threads);
  }

  return 0;
}
