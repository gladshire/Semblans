#include "paando.h"


void print_intro() {
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
  std::cout << std::left << std::setw(w.ws_col) << "                    P -ipeline for the       " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                    A -ssembly and           " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                    An -alysis of            " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                    D  -e novo                " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "       transcript-  O  -mics datasets         " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  ────────────────────────────────────────\n" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "A C++ package enabling the bulk retrieval," << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "assembly, and analysis of de novo transcriptomes" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "from multiple individuals\n" << std::endl;
}


void print_help() {
  std::cout << "  USAGE:\n" << std::endl;
  std::cout << "    paando [FUNCTION] [--config/-cfg] configuration_file_path "
            << "[--threads/-t] num_threads [--ram/-r] ram_to_dedicate\n" << std::endl;
}

int main(int argc, char * argv[]) {
  print_intro();
  // Command structure:
  //   ./paando preprocess/assemble/postprocess/{empty} config.ini threads ram_gb
  std::string threadStr;
  std::string ramStr;
  int numThreads;
  int ram;
  std::string pathConfig;
  std::string command;
  bool multAssembly;

  if (argc == 1 || 
      (argc == 2 && 
       (argv[1] == "--help" ||
        argv[1] == "-h"))) {
    print_help();
    exit(0);
  }
  else {
    // Parse through flags in commands
    threadStr = "";
    ramStr = "";
    numThreads = 1;
    ram = 1;
    std::string pathConfig = "";
    std::string command = "";
    multAssembly = false;
    
    for (int i = 0; i < argc; i++) {
      // Check for paando command (preprocess/assemble/postprocess)
      // If none is given, will perform all
      if (strcmp("preprocess", argv[i]) == 0 ||
          strcmp("Preprocess", argv[i]) == 0 ||
          strcmp("pre", argv[i]) == 0 ||
          strcmp("Pre", argv[i]) == 0 ||
          strcmp("pr", argv[i]) == 0 ||
          strcmp("Pr", argv[i]) == 0) {
        command = "preprocess";
      }
      else if (strcmp("assemble", argv[i]) == 0 ||
               strcmp("Assemble", argv[i]) == 0 ||
               strcmp("ass", argv[i]) == 0 ||
               strcmp("Ass", argv[i]) == 0 ||
               strcmp("as", argv[i]) == 0 ||
               strcmp("As", argv[i]) == 0 ||
               strcmp("a", argv[i]) == 0 ||
               strcmp("A", argv[i]) == 0) {
        command = "assemble";
      }
      else if (strcmp("postprocess", argv[i]) == 0 ||
               strcmp("Postprocess", argv[i]) == 0 ||
               strcmp("post", argv[i]) == 0 ||
               strcmp("Post", argv[i]) == 0 ||
               strcmp("pos", argv[i]) == 0 ||
               strcmp("Pos", argv[i]) == 0 ||
               strcmp("po", argv[i]) == 0 ||
               strcmp("Po", argv[i]) == 0) {
        command = "postprocess";
      }
      else if (strcmp("all", argv[i]) == 0 ||
               strcmp("All", argv[i]) == 0) {
        command = "all";
      }

      // Check for config file path
      else if (strcmp("--config", argv[i]) == 0 ||
               strcmp("-cfg", argv[i]) == 0) {
        if (strcmp(argv[i + 1] + strlen(argv[i + 1]) - 4, ".ini") == 0 ||
            strcmp(argv[i + 1] + strlen(argv[i + 1]) - 4, ".INI") == 0) {
          pathConfig = argv[i + 1];
          fs::path pathConfigFile(pathConfig.c_str());
          if (!fs::exists(pathConfigFile)) {
            // ERROR: Config file not found!
            exit(1);
          }
        }
        else {
          // ERROR: No config file specified!
          exit(1);
        }
      }

      // Check for number of threads flag
      else if (strcmp("--threads", argv[i]) == 0 ||
               strcmp("--Threads", argv[i]) == 0 ||
               strcmp("-t", argv[i]) == 0 ||
               strcmp("-T", argv[i]) == 0) {
        if (i != argc - 1) {
          threadStr = argv[i + 1];
          for (int j = 0; j < strlen(argv[i + 1]); j++) {
            if (!isdigit(threadStr[j])) {
              // ERROR: Bad thread argument
              exit(1);
            }
          }
          numThreads = atoi(threadStr.c_str());
        }
        else {
          // No thread count given, using 1 thread
        }
      }

      // Check for ammount of ram to be dedicated
      else if (strcmp("--ram", argv[i]) == 0 ||
               strcmp("--Ram", argv[i]) == 0 ||
               strcmp("-r", argv[i]) == 0 ||
               strcmp("-R", argv[i]) == 0) {
        if (i != argc - 1) {
          ramStr = argv[i + 1];
          for (int j = 0; j < strlen(argv[i + 1]); j++) {
            if (!isdigit(ramStr[j])) {
              // ERROR: Bad memory argument
              exit(1);
            }
          }
          ram = atoi(ramStr.c_str());
        }
        else {
          // No RAM ammount given, using 1 GB
        }
      }

      // Check if user specified a multi-SRA assembly
      else if (strcmp("--multi", argv[i]) == 0 ||
               strcmp("--Multi", argv[i]) == 0 ||
               strcmp("-m", argv[i]) == 0 ||
               strcmp("-M", argv[i]) == 0) {
        multAssembly = true;
      }
    }
    if (command == "") {
      command = "all";
    }
    if (pathConfig == "") {
      // ERROR: No config file specified
      exit(1);
    }
    // Case 1: preprocess
    std::string preCmd = PAANDO_DIR + "preprocess " + pathConfig + " " +
                         std::to_string(numThreads) + " " +
                         std::to_string(ram);
    std::string assCmd = PAANDO_DIR + "assemble " + pathConfig + " " +
                         std::to_string(numThreads) + " " +
                         std::to_string(ram);
    if (multAssembly) {
      assCmd += " --mult";
    }
    std::string postCmd = PAANDO_DIR + "postprocess " + pathConfig + " " +
                          std::to_string(numThreads) + " " +
                          std::to_string(ram);

    std::cout << preCmd << std::endl;
    // Case 1: preprocess
    if (command == "preprocess") {
      std::cout << "Performing preprocessing only ..." << std::endl;
      system(preCmd.c_str());
      exit(0);
    }
    // Case 2: assemble
    // TODO: Account for if there are no files for assembly
    if (command == "assemble") {
      std::cout << "Performing assembly only ..." << std::endl;
      system(assCmd.c_str());
      exit(0);
    }
    // Case 3: postprocess
    // TODO: Account for if there are no files for postprocess
    if (command == "postprocess") {
      std::cout << "Performing postprocess only ..." << std::endl;
      system(postCmd.c_str());
      exit(0);
    }
    // Case 4: all three
    if (command == "all") {
      system(preCmd.c_str());
      system(assCmd.c_str());
      system(postCmd.c_str());
      exit(0);
    }
  }
}
