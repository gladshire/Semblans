#include "postprocess.h"

// TODO:
//   Allow postprocessing of local transcripts
//   Otherwise retrieve from Trinity folder of project

std::vector<transcript> get_transcript() {
  std::vector<transcript> transcripts;
  fs::
  return transcripts;
}

void print_help() {

}

int main(int argc, char * argv[]) {
  if (argc > 1) {
    // Retrieve SRA objects, convert to transcripts
    INI_MAP cfgIni = make_ini_map(argv[1]);
    std::vector<SRA> sras;
    sras = get_sras(cfgIni);
    std::vector<transcript> transcripts = get_transcript(sras);
    // Get number of threads
    std::string threads = argv[2];
    // Get RAM in GB
    std::string ram_gb = argv[3];
    // Make blast db
    makeBlastDb("../uniprot-download_true_format_fasta_query__28_28proteome_3AUP00000080-2022.11.14-19.12.00.85.fasta", "../");
    // Run BlastX
     
  }
  else {
    print_help();
  }
}
