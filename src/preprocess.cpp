#include "preprocess.h"

void retrieve_sra_data(std::vector<SRA> sras, std::string threads) {
  std::cout << "Retrieving SRA runs for:\n" << std::endl;
  summarize_all_sras(sras);

  prefetch_sra(sras);
  fasterq_sra(sras, threads);
}


void print_help() {
  winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  std::cout << std::left << std::setw(w.ws_col) << "\n  ┌──────────────────────────────────────┐" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  │  _____                      _        │" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  │ |  __ \\                    | |       │" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  │ | |__) |_ _  __ _ _ __   __| | ___   │" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  │ |  ___/ _` |/ _` | '_ \\ / _` |/ _ \\  │ " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  │ | |  | (_| | (_| | | | | (_| | (_) | │" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  │ |_|   \\__,_|\\__,_|_| |_|\\__,_|\\___/  │" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  └──────────────────────────────────────┘\n" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                # #### ####" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "             ### \\/#|### |/####" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "           ##\\/#/ \\||/##/_/##/_#" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "         ###  \\/###|/ \\/ # ###" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "       ##_\\_#\\_\\## | #/###_/_####" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "      ## #### # \\ #| /  #### ##/##" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "       __#_--###`  |{,###---###-~" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                   \\-`/" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                    }P} -ipeline for the       " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                    /A{ -ssembly and           " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                   {AN} -alysis of            " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                   {D}  -e novo                " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "      transcript-  {O{  -mics datasets         " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  ────────────────────────────────────────\n" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "A C++ package enabling the bulk retrieval," << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "assembly, and analysis of de novo transcriptomes" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "from multiple individuals\n" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "COMMAND STRUCTURE" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "preprocess PATH/TO/CONFIG.INI num_threads RAM_GB" << std::endl;
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
    std::vector<SRA> sras;
    std::string localDataBoolStr = cfgIni["General"]["use_local_data"];
    bool localDataBool;
    for (int i = 0; i < localDataBoolStr.length(); i++) {
      localDataBoolStr[i] = std::tolower(localDataBoolStr[i]);
    }
    localDataBool = (localDataBoolStr == "true") ? true : false;
    // TODO: Add error message for value that isn't 'true' or 'false'
    if (localDataBool) {
      std::string localDataDirStr = cfgIni["General"]["local_data_directory"];
      fs::path localDataDir(localDataDirStr.c_str());
      fs::directory_iterator localDirIter{localDataDir};
      std::vector<fs::path> filesSkip;
      while (localDirIter != fs::directory_iterator{}) {
        fs::path sraFile = localDirIter->path();
        if (std::find(filesSkip.begin(), filesSkip.end(), sraFile) != filesSkip.end()) {
          localDirIter++;
          continue;
        }
        std::string sraFileName = sraFile.filename().c_str();
        std::string sraFileNameP = sraFileName;
        std::string sraFile1;
        std::string sraFile2;
        if (std::string(sraFile.stem().c_str()).back() == '1' ||
            std::string(sraFile.stem().c_str()).back() == '2') {
          if (std::string(sraFile.stem().c_str()).back() == '1') {
            sraFileNameP.replace(sraFileName.length() - 7, 1, "2");
            sraFile1 = sraFileName;
            sraFile2 = sraFileNameP;
          }
          else if (std::string(sraFile.stem().c_str()).back() == '2') {
            sraFileNameP.replace(sraFileName.length() - 7, 1, "1");
            sraFile1 = sraFileNameP;
            sraFile2 = sraFileName;
          }
          else {
            // File not paired-end
            //   Construct UNPAIRED SRA object with filename
            //   Add filename to filesskip
            SRA currSra(sraFileName, "", cfgIni);
            sras.push_back(currSra);
            filesSkip.push_back(sraFile);
            continue;
          }
          fs::path sraFileP = sraFile.parent_path() / fs::path(sraFileNameP);
          if (fs::exists(sraFileP)) {
            SRA currSra(sraFile1, sraFile2, cfgIni);
            sras.push_back(currSra);
            filesSkip.push_back(sraFile);
            filesSkip.push_back(sraFileP);
          }
          else {
            SRA currSra(sraFileName, "", cfgIni);
            sras.push_back(currSra);
            filesSkip.push_back(sraFile);
          }
        }
        else {
          // File not paired-end
          SRA currSra(sraFileName, "", cfgIni);
          sras.push_back(currSra);
          filesSkip.push_back(sraFile);
        }
        localDirIter++;
      }
    }
    else {
      // Not using local files, retrieve from NCBI instead
      sras = get_sras(cfgIni);
    }
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
