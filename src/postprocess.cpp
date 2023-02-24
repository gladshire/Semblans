#include "postprocess.h"

// TODO:
//   Allow postprocessing of local transcripts
//   Otherwise retrieve from Trinity folder of project

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
    make_proj_space(cfgIni);
    std::string outputDir = cfgIni["General"]["output_directory"];
    std::string projDir = outputDir + cfgIni["General"]["project_name"] + "/";
    std::string refProt = cfgIni["General"]["reference_proteome_path"];

    std::vector<SRA> sras;
    std::vector<std::string> localDataFiles;
    sras = get_sras(cfgIni);
    for (auto fqFileName : cfgIni.at("Local files")) {
      localDataFiles.push_back(fqFileName.first);
    }
    std::pair<std::string, std::string> sraRunsLocal;
    size_t pos;
    for (auto sraRun : localDataFiles) {
      sraRunsLocal.first = "";
      sraRunsLocal.second = "";
      pos = sraRun.find(" ");
      sraRunsLocal.first = sraRun.substr(0, pos);
      if (pos != std::string::npos) {
        sraRun.erase(0, pos + 1);
        pos = sraRun.find(" ");
        sraRunsLocal.second = sraRun.substr(0, pos);
      }
      if (fs::exists(cfgIni["General"]["local_data_directory"] + sraRunsLocal.first) &&
          fs::exists(cfgIni["General"]["local_data_directory"] + sraRunsLocal.second)) {
        sras.push_back(SRA(sraRunsLocal.first, sraRunsLocal.second, cfgIni));
      }
      else {
        if (sraRunsLocal.first != "" &&
            !fs::exists(cfgIni["General"]["local_data_directory"] + sraRunsLocal.first)) {
          std::cout << "ERROR: Local run not found: \"" << sraRunsLocal.first << "\""
                    << std::endl;
        }
        if (sraRunsLocal.second != "" &&
            !fs::exists(cfgIni["General"]["local_data_directory"] + sraRunsLocal.second)) {
          std::cout << "ERROR: Local run not found: \"" << sraRunsLocal.second << "\""
                    << std::endl;
        }
      }
    }
    transcript trans = get_transcript_mult(sras[0]);
    // Get number of threads
    std::string threads = argv[2];
    // Get RAM in GB
    std::string ram_gb = argv[3];
    // Make blast db
    makeDb(refProt, projDir + stepDirs[8]);
    // Run BlastX
    std::string blastDbName = refProt.substr(refProt.find_last_of("/"),
                                             refProt.find_last_of(".") -
                                             refProt.find_last_of("/"));
    blastxDiam(trans, projDir + stepDirs[8] + blastDbName, threads,
               projDir + stepDirs[8]);
    // Detect and remove chimeric transcripts
    detect_chimera(trans, std::string(trans.get_trans_path_blastx().c_str()),
                   projDir + stepDirs[8]);
    removeChimera(trans, std::string(trans.get_trans_path_cinfo().c_str()),
                         std::string(trans.get_trans_path_ccut().c_str()),
                         (uintmax_t)(stoi(ram_gb) * 1000000000),
                         std::string(trans.get_trans_path_chimera().parent_path().c_str()));
    // Perform salmon index of transcripts
    salmon_index(trans, threads);
    // Perform salmon quant of transcripts / reads
    salmon_quant(trans, sras, threads);
    // Perform corset run to cluster transcripts
    corset_eq_classes(trans, std::string((trans.get_trans_path_quant() /
                                         fs::path("aux_info") /
                                         fs::path("eq_classes.txt.gz")).c_str()),
                             std::string(projDir + stepDirs[9]));
    // Filter corset output
    filterCorset(trans, std::string(trans.get_trans_path_clust().c_str()),
                 (uintmax_t)(stoi(ram_gb) * 1000000000),
                 std::string(projDir + stepDirs[9]));
    // Run transdecoder
    run_transdecoder(trans, threads, (uintmax_t)(stoi(ram_gb) * 1000000000),
                     projDir + stepDirs[8] + blastDbName, projDir + stepDirs[10]);
  }
  else {
    print_help();
  }
}
