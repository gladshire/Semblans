#include "preprocess.h"

std::vector<std::string> stepDirs = {"00-Raw_reads/", "01-Quality_analysis_1/",
                                     "02-Error_correction/", "03-Trimming/",
                                     "04-Filter_foreign/", "05-Quality_analysis_2/",
                                     "06-Filter_overrepresented/"};


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
  std::cout << "\n" << "preprocess PATH/TO/CONFIG.INI num_threads RAM_GB" << std::endl;
}


int main(int argc, char * argv[]) {
  if (argc != 4) {
    print_help();
    return 0;
  }
  else {
    INI_MAP cfgIni = make_ini_map(argv[1]);
    std::string threads = argv[2];
    std::string ram_gb = argv[3];
    std::vector<SRA> sras = get_sras(cfgIni);
    std::string fastqc_dir_1(sras[0].get_fastqc_dir_1().first.parent_path().parent_path().c_str());
    std::string fastqc_dir_2(sras[0].get_fastqc_dir_2().first.parent_path().parent_path().c_str());
    std::vector<std::string> kraken2Dbs = get_kraken2_dbs(cfgIni);
    std::string kraken2_conf = get_kraken2_conf(cfgIni);
    std::pair<std::vector<std::string>, std::vector<std::string>> overrepSeqs;
    make_proj_space(cfgIni);
    retrieve_sra_data(sras, threads);
    run_fastqc_bulk(sras, threads, fastqc_dir_1);
    run_rcorr(sras, threads);
    rem_unfix_bulk(sras, ram_gb);
    run_trimmomatic(sras, threads);
    run_kraken2_dbs(sras, threads, kraken2Dbs, kraken2_conf);
    run_fastqc_bulk(sras, threads, fastqc_dir_2);
    rem_overrep_bulk(sras, ram_gb);
  }

  return 0;
}
