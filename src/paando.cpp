#include "paando.h"


void print_intro(std::string logFile) {
  winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  logOutput("\n  ┌──────────────────────────────────────┐", logFile);
  logOutput("  │  _____                      _        │", logFile);
  logOutput("  │ |  __ \\                    | |       │", logFile);
  logOutput("  │ | |__) |_ _  __ _ _ __   __| | ___   │", logFile);
  logOutput("  │ |  ___/ _` |/ _` | '_ \\ / _` |/ _ \\  │ ", logFile);
  logOutput("  │ | |  | (_| | (_| | | | | (_| | (_) | │", logFile);
  logOutput("  │ |_|   \\__,_|\\__,_|_| |_|\\__,_|\\___/  │", logFile);
  logOutput("  └──────────────────────────────────────┘\n", logFile);
  logOutput("                    P  -ipeline for the       ", logFile);
  logOutput("                    A  -ssembly and           ", logFile);
  logOutput("                    An -alysis of            ", logFile);
  logOutput("                    D  -e novo                ", logFile);
  logOutput("       transcript-  O  -mics datasets         ", logFile);
  logOutput("  ────────────────────────────────────────\n", logFile);
  logOutput("A C++ package enabling the bulk retrieval,", logFile);
  logOutput("assembly, and analysis of de novo transcriptomes", logFile);
  logOutput("from multiple individuals\n", logFile);
  logOutput("  ────────────────────────────────────────\n", logFile);
}


void print_help() {
  std::cout << "USAGE:\n" << std::endl;
  std::cout << "  paando [--help/-h] [COMMAND] [--config/-cfg]\n"
            << "         [--threads/-t] [--ram/-r] [--multi/-m]\n"
            << "         [--retain/-f] [--verbose/-v]\n" << std::endl;
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
  std::cout << "  -f,   --retain     Prevents deletion of intermediate files in pipeline" << std::endl;
  std::cout << "  -v,   --verbose    Prints all output from Paando and sub-programs" << std::endl;
  std::cout << "  -h,   --help       Displays this help screen" << std::endl;

}

int main(int argc, char * argv[]) {
  //print_intro();
  std::string threadStr;
  std::string ramStr;
  int numThreads;
  int ram;
  std::string pathConfig;
  std::string command;
  bool multAssembly;
  bool useLocalData;
  bool retainInterFiles;
  bool verboseOutput;

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
    verboseOutput = false;

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

      else if (strcmp("--retain", argv[i]) == 0 ||
               strcmp("--Retain", argv[i]) == 0 ||
               strcmp("-f", argv[i]) == 0 ||
               strcmp("-F", argv[i]) == 0) {
        retainInterFiles = true;
      }
      else if (strcmp("--verbose", argv[i]) == 0 ||
               strcmp("--Verbose", argv[i]) == 0 ||
               strcmp("-v", argv[i]) == 0 ||
               strcmp("-V", argv[i]) == 0) {
        verboseOutput = true;
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

    INI_MAP cfgIni = make_ini_map(pathConfig.c_str());
    std::string logFilePath = cfgIni["General"]["log_file"];
    std::ofstream logFile(logFilePath, std::ios_base::trunc);
 
    std::string verbose;
    std::string retain;
    if (verboseOutput) {
      verbose = " true";
    }
    else {
      verbose = " false";
    }
    if (retainInterFiles) {
      retain = " true";
    }
    else {
      retain = " false";
    }
    std::string preCmd = PAANDO_DIR + "preprocess " + pathConfig + " " +
                         std::to_string(numThreads) + " " +
                         std::to_string(ram);
    preCmd += retain;
    preCmd += verbose;
    std::string assCmd = PAANDO_DIR + "assemble " + pathConfig + " " +
                         std::to_string(numThreads) + " " +
                         std::to_string(ram);
    if (multAssembly) {
      assCmd += " true";
    }
    else {
      assCmd += " false";
    }
    assCmd += retain;
    assCmd += verbose;
    std::string postCmd = PAANDO_DIR + "postprocess " + pathConfig + " " +
                          std::to_string(numThreads) + " " +
                          std::to_string(ram);
    postCmd += retain;
    postCmd += verbose;
    int result;

    print_intro(logFilePath);

    // Case 1: preprocess
    if (command == "preprocess") {
      logOutput("Performing preprocessing only", logFilePath);
      result = system(preCmd.c_str());
      if (WIFSIGNALED(result)) {
        exit(1);
      }
      exit(0);
    }
    // Case 2: assemble
    if (command == "assemble") {
      logOutput("Performing assembly only", logFilePath);
      result = system(assCmd.c_str());
      if (WIFSIGNALED(result)) {
        exit(1);
      }
      exit(0);
    }
    // Case 3: postprocess
    if (command == "postprocess") {
      logOutput("Performing postprocess only", logFilePath);
      result = system(postCmd.c_str());
      if (WIFSIGNALED(result)) {
        exit(1);
      }
      exit(0);
    }
    // Case 4: all three
    if (command == "all") {
      logOutput("Performing entire assembly", logFilePath);
      result = system(preCmd.c_str());
      if (WIFSIGNALED(result)) {
        exit(1);
      }
      logOutput("", logFilePath);
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
