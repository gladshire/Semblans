#include "postprocess.h"

// TODO:
//   Allow postprocessing of local transcripts
//   Otherwise retrieve from Trinity folder of project

std::vector<transcript> get_transcript(std::vector<SRA> sras) {
  std::vector<transcript> transcripts;
  for (auto & sra : sras) {
    transcript currTrans(sra);
    transcripts.push_back(currTrans);
  }
  return transcripts;
}

transcript get_transcript_mult(SRA sra) {
  transcript trans(sra);
  return trans;
}

void print_help() {

}

int main(int argc, char * argv[]) {
  if (argc > 1) {
    // Retrieve SRA objects, convert to transcripts
    INI_MAP cfgIni = make_ini_map(argv[1]);
    std::vector<SRA> sras;
    sras = get_sras(cfgIni);
    transcript trans = get_transcript_mult(sras[0]);
    // Get number of threads
    std::string threads = argv[2];
    // Get RAM in GB
    std::string ram_gb = argv[3];
    // Make blast db
    makeBlastDb("../uniprot-download_true_format_fasta_query__28_28proteome_3AUP00000080-2022.11.14-19.12.00.85.fasta", "../");
    // Run BlastX
    blastx(trans, "../uniprot-download_true_format_fasta_query__28_28proteome_3AUP00000080-2022.11.14-19.12.00.85", threads, "../");
    detect_chimera(trans, "../7227_Drosophila_melanogaster.Trinity.blastx", "../");
    removeChimera(trans, "../7227_Drosophila_melanogaster.Trinity.info",
                         "../7227_Drosophila_melanogaster.Trinity.cut",
                         (uintmax_t)(stoi(ram_gb) * 1000000000),
                         ".");
  }
  else {
    print_help();
  }
}
