#include "assemble.h"

// TODO:
//   Allow assembly of local short reads
//   Otherwise retrieve from overrepresented folder

std::vector<transcript> get_transcript(std::vector<SRA> sras) {
  std::vector<transcript> transcripts;
  for (auto &sra : sras) {
    transcript currTrans(sra);
    transcripts.push_back(currTrans);
  }
  return transcripts;
}

void print_help() {
  std::cout << "\n" << "NAME_OF_PROGRAM" << " - "
            << "A tool for bulk assemblies of de novo transcriptome data" << std::endl;
  std::cout << "\n" << "COMMAND STRUCTURE" << std::endl;
  std::cout << "\n" << "assemble PATH/TO/CONFIG.INI num_threads RAM_GB (--mult)" << std::endl;
}

int main(int argc, char * argv[]) {
  if (argc > 1) {
    std::vector<SRA> sras;
    // Retrieve SRA objects for trinity runs
    INI_MAP cfgIni = make_ini_map(argv[1]);
    sras = get_sras(cfgIni);
    // Get number of threads
    std::string threads = argv[2];
    // Get RAM in GB
    std::string ram_gb = argv[3];
    // Get boolean for multiple sra processing
    std::string mult_sra_str;
    bool mult_sra = false;
    if (argc == 5) {
      mult_sra_str = argv[4];
      if (mult_sra_str == "--mult") {
        mult_sra = true;
      }
    }
    // Perform assembly with Trinity
    std::vector<transcript> transcriptsSra = run_trinity_bulk(sras, threads, ram_gb, mult_sra);
  }
  else {
    print_help();
  }
  return 0;
}
