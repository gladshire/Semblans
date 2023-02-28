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
  std::cout << std::left << std::setw(w.ws_col) << "                    P  -ipeline for the       " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                    A  -ssembly and           " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                    An -alysis of            " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "                    D  -e novo                " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "       transcript-  O  -mics datasets         " << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  ────────────────────────────────────────\n" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "A C++ package enabling the bulk retrieval," << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "assembly, and analysis of de novo transcriptomes" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "from multiple individuals\n" << std::endl;
  std::cout << std::left << std::setw(w.ws_col) << "  ────────────────────────────────────────\n" << std::endl;

}


void print_help() {
  std::cout << "USAGE:\n" << std::endl;
  std::cout << "  paando [--help/-h] [COMMAND] [--config/-cfg]\n"
            << "         [--threads/-t] [--ram/-r] [--multi/-m]\n" << std::endl;
  std::cout << "ARGUMENTS:\n" << std::endl;
  std::cout << "  [COMMAND]" << std::endl;
  std::cout << "    preprocess       Performs pre-assembly steps only" << std::endl;
  std::cout << "    assemble         Performs de novo assembly step only" << std::endl;
  std::cout << "    postprocess      Performs post-assembly steps only" << std::endl;
  std::cout << "    all (default)    Performs all steps in pipeline\n" << std::endl;
  std::cout << "  -cfg, --config     Specifies path to configuration file (REQUIRED)" << std::endl;
  std::cout << "  -t,   --threads    Specifies number of threads/CPU cores to employ" << std::endl;
  std::cout << "  -r,   --ram        Specifies ammount of memory/RAM (GB) to dedicate" << std::endl;
  std::cout << "  -m,   --multi      Perform assembly from multiple SRA runs" << std::endl;
  std::cout << "  -h,   --help       Displays this help screen" << std::endl;

}

int main(int argc, char * argv[]) {
  print_intro();
  std::string threadStr;
  std::string ramStr;
  int numThreads;
  int ram;
  std::string pathConfig;
  std::string command;
  bool multAssembly;
  bool useLocalData;
  bool retainInterFiles;

  if (argc == 1 || 
      (argc == 2 && 
       (strcmp("--help", argv[1]) == 0 ||
        strcmp("-h", argv[1]) == 0))) {
    print_help();
    exit(0);
  }
  if (argc == 2 && (strcmp("--version", argv[1]) == 0 ||
                    strcmp("-v", argv[1]) == 0)) {
    std::cout << "Paando version: 0.0.1" << std::endl;
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
    retainInterFiles = false;    

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
            std::cout << "ERROR: Config file: " << pathConfig << " not found\n" << std::endl;
            exit(1);
          }
        }
        else {
          // ERROR: No config file specified!
          std::cout << "ERROR: Must specify config file" << std::endl;
          std::cout << "  (example: --config path/to/config.ini)\n" << std::endl;
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
              std::cout << "ERROR: Invalid argument given for --threads/-t" << std::endl;
              std::cout << "  (example: --threads 8)\n" << std::endl;
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
              std::cout << "ERROR: Invalid argument given for --ram/-r" << std::endl;
              std::cout << "  (example: --ram 10)\n" << std::endl;
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
               strcmp("--mult", argv[i]) == 0 ||
               strcmp("--Mult", argv[i]) == 0 ||
               strcmp("-m", argv[i]) == 0 ||
               strcmp("-M", argv[i]) == 0) {
        multAssembly = true;
      }

      else if (strcmp("--retain-intermediates", argv[i]) == 0 ||
               strcmp("--Retain-intermediates", argv[i]) == 0 ||
               strcmp("--retain-Intermediates", argv[i]) == 0 ||
               strcmp("--Retain-Intermediates", argv[i]) == 0 ||
               strcmp("-i", argv[i]) == 0 ||
               strcmp("-I", argv[i]) == 0) {
        retainInterFiles = true;
      }
    }
    if (command == "") {
      command = "all";
    }
    if (pathConfig == "") {
      // ERROR: No config file specified
      std::cout << "ERROR: Must specify config file" << std::endl;
      std::cout << "  (example: --config path/to/config.ini)\n" << std::endl;
      exit(1);
    }
    std::string preCmd = PAANDO_DIR + "preprocess " + pathConfig + " " +
                         std::to_string(numThreads) + " " +
                         std::to_string(ram);
    if (retainInterFiles) {
      preCmd += " true";
    }
    else {
      preCmd += " false";
    }
    std::string assCmd = PAANDO_DIR + "assemble " + pathConfig + " " +
                         std::to_string(numThreads) + " " +
                         std::to_string(ram);
    if (multAssembly) {
      assCmd += " --mult";
    }
    std::string postCmd = PAANDO_DIR + "postprocess " + pathConfig + " " +
                          std::to_string(numThreads) + " " +
                          std::to_string(ram);
    int result;
    // Case 1: preprocess
    if (command == "preprocess") {
      std::cout << "Performing preprocessing only ..." << std::endl;
      result = system(preCmd.c_str());
      if (WIFSIGNALED(result)) {
        exit(1);
      }
      exit(0);
    }
    // Case 2: assemble
    // TODO: Account for if there are no files for assembly
    if (command == "assemble") {
      std::cout << "Performing assembly only ..." << std::endl;
      result = system(assCmd.c_str());
      if (WIFSIGNALED(result)) {
        exit(1);
      }
      exit(0);
    }
    // Case 3: postprocess
    // TODO: Account for if there are no files for postprocess
    if (command == "postprocess") {
      std::cout << "Performing postprocess only ..." << std::endl;
      result = system(postCmd.c_str());
      if (WIFSIGNALED(result)) {
        exit(1);
      }
      exit(0);
    }
    // Case 4: all three
    if (command == "all") {
      result = system(preCmd.c_str());
      if (WIFSIGNALED(result)) {
        exit(1);
      }
      result = system(assCmd.c_str());
      if (WIFSIGNALED(result)) {
        exit(1);
      }
      result = system(postCmd.c_str());
      if (WIFSIGNALED(result)) {
        exit(1);
      }
      exit(0);
    }
  }
}
