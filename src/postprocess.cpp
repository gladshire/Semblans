#include "postprocess.h"

// TODO:
//   Implemente bulk postprocess of multiple transcripts
//   Currently, only performs for single transcript
//   Retrieve from Trinity folder of project

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

bool stringToBool(std::string boolStr) {
  bool boolConv;
  for (int i = 0; i < boolStr.length(); i++) {
    boolStr[i] = std::tolower(boolStr[i]);
  }
  boolConv = (boolStr == "true") ? true : false;
  return boolConv;
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
    std::string logFilePath = cfgIni["General"]["log_file"];
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
    if (sras.empty()) {
      std::cout << "ERROR: No SRA runs specified. Please check config file" << std::endl;
    }
    // TODO: Pipeline only performs postprocess on one transcript
    //   Make it instantiate transcript object for each file in trinity folder
    transcript trans = get_transcript_mult(sras[0]);
    std::vector<transcript> transVec;
    transVec.push_back(trans);



    // Get number of threads
    std::string threads = argv[2];
    // Get RAM in GB
    std::string ram_gb = argv[3];
    uintmax_t ram_b = (uintmax_t)(stoi(ram_gb)) * 1000000000;
    // Determine whether to print output of programs
    bool dispOutput = stringToBool(argv[4]);

    // Summarize program execution parameters
    logOutput("Paando Postprocess started with following parameters:", logFilePath);
    logOutput("  Config file:     " + std::string(argv[1]), logFilePath);
    logOutput("  Threads (Cores): " + threads, logFilePath);
    logOutput("  Memory (GB):     " + ram_gb, logFilePath);
    logOutput("  Reference Prot:  " + refProt, logFilePath);
    logOutput("  SRA runs:\n", logFilePath);
    summarize_all_sras(sras, logFilePath, 6);


    // BlastX alignment of transcript to reference proteome
    std::string currTransInDiam;
    std::string currBlastDbName;

    std::string blastDbName = refProt.substr(refProt.find_last_of("/"),
                                             refProt.find_last_of(".") -
                                             refProt.find_last_of("/"));
    for (auto trans : transVec) {
      currTransInDiam = trans.get_trans_path_trinity().c_str();

      makeDb(refProt, projDir + stepDirs[8], dispOutput, logFilePath);
      // Run BlastX
      currBlastDbName = refProt.substr(refProt.find_last_of("/"),
                                       refProt.find_last_of(".") -
                                       refProt.find_last_of("/"));
      blastxDiam(currTransInDiam, projDir + stepDirs[8] + currBlastDbName, threads,
                 projDir + stepDirs[8], dispOutput, logFilePath);
    }

    // Detect and remove chimeric transcripts
    std::string currTransIn;
    std::string currTransOut;
    std::string currBlastx;
    std::string currTransInfo;
    std::string currTransCut;
    std::string chimOutDir = projDir + stepDirs[8];
    for (auto trans : transVec) {
      currTransIn = trans.get_trans_path_trinity().c_str();
      currTransOut = trans.get_trans_path_chimera().c_str();

      currBlastx = trans.get_trans_path_blastx().c_str();
      currTransInfo = trans.get_trans_path_cinfo().c_str();
      currTransCut = trans.get_trans_path_ccut().c_str();

      detect_chimera(currBlastx, chimOutDir);
      removeChimera(currTransIn, currTransOut, currTransInfo, currTransCut, ram_gb,
                    logFilePath);
    }
    // Perform salmon index of transcripts
    salmon_index(trans, threads, dispOutput, logFilePath);
    // Perform salmon quant of transcripts / reads
    salmon_quant(trans, sras, threads, dispOutput, logFilePath);
    // Perform corset run to cluster transcripts
    corset_eq_classes(trans, std::string((trans.get_trans_path_quant() /
                                         fs::path("aux_info") /
                                         fs::path("eq_classes.txt.gz")).c_str()),
                             std::string(projDir + stepDirs[9]), dispOutput, logFilePath);
    // Filter corset output
    filterCorset(trans, std::string(trans.get_trans_path_clust().c_str()),
                 ram_b,
                 std::string(projDir + stepDirs[9]),
                 logFilePath);
    // Run transdecoder
    run_transdecoder(trans, threads, ram_b,
                     projDir + stepDirs[8] + blastDbName, projDir + stepDirs[10],
                     dispOutput, logFilePath);
  }
  else {
    print_help();
  }
}
