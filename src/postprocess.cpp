#include "postprocess.h"

// TODO:
//   Allow postprocessing of local transcripts
//   Otherwise retrieve from Trinity folder of project

// TODO:
//   Blast wrapping   (X)
//   Chimera filter   (X)
//   Salmon wrapping  (X)
//   Corset wrapping  (X)
//   Corset filtering ( )
//   Transdecoder wrp ( )

extern std::vector<std::string> stepDirs;

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
    std::string ref_prot_path = cfgIni["General"]["reference_proteome_path"];
    makeBlastDb(ref_prot_path, cfgIni["General"]["output_directory"] + cfgIni["General"]["project_name"] + "/" + stepDirs[8]);
    // Run BlastX
    std::string blastDbName = ref_prot_path.substr(ref_prot_path.find_last_of("/"),
                                                   ref_prot_path.find_last_of(".") - ref_prot_path.find_last_of("/"));
    blastx(trans, cfgIni["General"]["output_directory"] + cfgIni["General"]["project_name"] + "/" + stepDirs[8] + blastDbName, threads, "../");
    // TODO: Make the below commands' directories dynamic
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
