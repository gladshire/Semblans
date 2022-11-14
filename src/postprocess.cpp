#include "postprocess.h"

std::vector<transcript> get_transcript(std::vector<SRA> sras) {
  std::vector<transcript> transcripts;
  for (auto &sra : sras) {
    transcript currTrans(sra);
    transcripts.push_back(currTrans);
  }
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
  }
  else {
    print_help();
  }
}
