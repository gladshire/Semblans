#include "semblans.h"

std::atomic<bool> procRunning(false);

// Print ASCII startup screen for Semblans
void print_intro(std::string logFile) {
  winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
  logOutput("\n", logFile);
  logOutput("\n", logFile);
  logOutput("                            |     |                     \n", logFile);
  logOutput("      __|   _ \\  __ `__ \\   __ \\  |   _` |  __ \\    __| \n", logFile);
  logOutput("    \\__ \\   __/  |   |   |  |   | |  (   |  |   | \\__ \\ \n", logFile);
  logOutput("    ____/ \\___| _|  _|  _| _.__/ _| \\__,_| _|  _| ____/ \n", logFile); 
  logOutput("\n\n\n", logFile);
  printBreakLine(logFile, 6, 47);
}

// Print base help message for Semblans in terminal
void print_help_base() {
  std::cout << "USAGE:\n" << std::endl;
  std::cout << "  semblans [--help/-h] [COMMAND] [--config/-cfg]\n"
            << "           [--threads/-t] [--ram/-r] [--retain/-f]\n"
            << "           [--verbose/-v]\n" << std::endl;
  std::cout << "ARGUMENTS:\n" << std::endl;
  std::cout << "  [COMMAND]" << std::endl;
  std::cout << "    preprocess           Performs pre-assembly steps only" << std::endl;
  std::cout << "    assemble             Performs de novo assembly step only" << std::endl;
  std::cout << "    postprocess          Performs post-assembly steps only" << std::endl;
  std::cout << "    all (default)        Performs all steps in pipeline\n" << std::endl;
  std::cout << "  [MODE]" << std::endl;
  std::cout << "  -1/-2, --left/--right  Specifies simple mode." << std::endl;
  std::cout << "                         Path to left/right read FASTQ(s) should follow the corresponding flags." << std::endl;
  std::cout << "                         To specify multiple pairs of files, separate file paths with commas." << std::endl;
  std::cout << "  -o, --out-directory   Specifies path to where Semblans should place output files." << std::endl;
  std::cout << "  -p, --prefix           Specifies the name of the resulting assembly Semblans creates." << std::endl;
  std::cout << "                         Directory must already exist" << std::endl;
  std::cout << "  -kdb, --kraken-db      Specifies which Kraken databases to use during removal of foreign reads." << std::endl;
  std::cout << "                         To specify multiple databases, separate their paths with commas." << std::endl;
  std::cout << "                         If omitted, Semblans will skip foreign read filtration" << std::endl;
  std::cout << "  -a, --assembly         Specifies an assembly for Semblans to process." << std::endl;
  std::cout << "                         Only used during postprocessing." << std::endl;
  std::cout << "  -rp, --ref-proteome    Specifies the reference proteome Semblans will use for" << std::endl;
  std::cout << "                         chimera detection and coding region prediction." << std::endl;
  std::cout << "                         Only used during postprocessing." << std::endl;
  std::cout << std::endl;
  std::cout << "  -cfg, --config         Specifies config file mode." << std::endl;
  std::cout << "                         Path to config file should follow this argument" << std::endl;
  std::cout << std::endl;
  std::cout << "  [FLAGS]" << std::endl;
  std::cout << "  -t,   --threads        Specifies number of threads/CPU cores to employ" << std::endl;
  std::cout << "  -r,   --ram            Specifies ammount of memory/RAM (GB) to dedicate" << std::endl;
  std::cout << "  -f,   --retain         Prevents deletion of intermediate files in pipeline" << std::endl;
  std::cout << "  -v,   --verbose        Prints all output from Semblans and sub-programs" << std::endl;
  std::cout << "  -h,   --help           Displays this help screen" << std::endl;

}


void print_help_preprocess() {
  std::cout << "USAGE:\n" << std::endl;
  std::cout << "  semblans preprocess [--help/-h] [--config/-cfg]\n"
            << "                      [--left/-l reads_1_left.fq,reads_2_left.fq,...]\n"
            << "                      [--right/-r reads_1_right.fq,reads_2_right.fq,...]\n"
            << "                      [--threads/-t] [--ram/-r] [--retain/-f]\n"
            << "                      [--verbose/-v]\n" << std::endl;
  std::cout << "ARGUMENTS:\n" << std::endl;
  std::cout << "  -1/-2, --left/--right  Specifies simple mode." << std::endl;
  std::cout << "                         Path to left/right read FASTQ(s) should follow the corresponding flags." << std::endl;
  std::cout << "                         To specify multiple pairs of files, separate file paths with commas" << std::endl;
  std::cout << "  -o,  --out-directory  Specifies path to where Semblans should place output files." << std::endl;
  std::cout << "  -p, --prefix           Specifies the name of the resulting assembly Semblans creates." << std::endl;
  std::cout << "                         Directory must already exist" << std::endl;
  std::cout << "  -kdb, --kraken-db      Specifies which Kraken databases to use during removal of foreign reads." << std::endl;
  std::cout << "                         To specify multiple databases, separate their paths with commas." << std::endl;
  std::cout << "  -cfg, --config         Specifies config file mode." << std::endl;
  std::cout << "                         Path to config file should follow this argument" << std::endl;
  std::cout << "  -t,   --threads        Specifies number of threads/CPU cores to employ" << std::endl;
  std::cout << "  -r,   --ram            Specifies ammount of memory/RAM (GB) to dedicate" << std::endl;
  std::cout << "  -f,   --retain         Prevents deletion of intermediate files in pipeline" << std::endl;
  std::cout << "  -v,   --verbose        Prints all output from Semblans and sub-programs" << std::endl;
  std::cout << "  -h,   --help           Displays this help screen" << std::endl;
}

void print_help_postprocess() {
  std::cout << "USAGE:\n" << std::endl;
  std::cout << "  semblans postprocess [--help/-h] [--config/-cfg]\n"
            << "                       [--assembly/-a reads_assembly.fa]\n"
            << "                       [--left/-1 reads_1_left.fq,reads_2_left.fq,...]\n"
            << "                       [--right/-2 reads_1_right.fq,reads_2_right.fq,...]\n"
            << "                       [--threads/-t] [--ram/-r] [--retain/-f]\n"
            << "                       [--verbose/-v]\n" << std::endl;
  std::cout << "ARGUMENTS:\n" << std::endl;
  std::cout << "  -1/-2, --left/--right  Specifies simple mode." << std::endl;
  std::cout << "                         Path to left/right read FASTQ(s) should follow the corresponding flags." << std::endl;
  std::cout << "                         To specify multiple pairs of files, separate file paths with commas" << std::endl;
  std::cout << "  -a,   --assembly       Specifies an assembly for Semblans to operate on." << std::endl;
  std::cout << "                         Only used during postprocessing." << std::endl;
  std::cout << "  -o,  --out-directory  Specifies path to where Semblans should place output files." << std::endl;
  std::cout << "  -p, --prefix           Specifies the name of the resulting assembly Semblans creates." << std::endl;
  std::cout << "                         Directory must already exist" << std::endl;
  std::cout << "  -cfg, --config         Specifies config file mode." << std::endl;
  std::cout << "                         Path to config file should follow this argument" << std::endl;
  std::cout << "  -t,   --threads        Specifies number of threads/CPU cores to employ" << std::endl;
  std::cout << "  -r,   --ram            Specifies ammount of memory/RAM (GB) to dedicate" << std::endl;
  std::cout << "  -f,   --retain         Prevents deletion of intermediate files in pipeline" << std::endl;
  std::cout << "  -v,   --verbose        Prints all output from Semblans and sub-programs" << std::endl;
  std::cout << "  -h,   --help           Displays this help screen" << std::endl;
 
}

// Utility function to obtain name of cleaned read files from raw file names
std::vector<std::string> makeCleanedNames(std::string commaSepReadFiles) {
  std::vector<std::string> cleanedReads;
  fs::path currFilePath;
  size_t commaInd;
  size_t currPos = 0;
  std::string currStr;
  do {
    commaInd = commaSepReadFiles.find(',', currPos);
    currStr = commaSepReadFiles.substr(currPos, commaInd - currPos);
    currFilePath = fs::path(currStr.c_str()).stem();
    currFilePath.replace_extension(".orep.filt.fq");
    cleanedReads.push_back(std::string(currFilePath.c_str()));
    currPos = commaInd + 1;
  } while (commaInd != std::string::npos);
  return cleanedReads;
}


void parseArgv(int argc, char * argv[], std::string & command,
               std::string & leftReads, std::string & rightReads,
               std::string & assembly, std::string & refProt,
               std::string & outDir, std::string & outPrefix,
               std::string & threadStr, std::string & ramStr,
               std::string & pathConfig, std::string & kraken2Dbs,
               bool & retainInterFiles, bool & verboseOutput) {

  bool argIsFlag;
  std::vector<int> nonFlagInd;
  nonFlagInd.push_back(0);
  for (int i = 0; i < argc; i++) {
    argIsFlag = true;
    // Check for semblans command (preprocess/assemble/postprocess)
    // If none is given, will perform all
    if (strcmp("preprocess", argv[i]) == 0 ||
        strcmp("Preprocess", argv[i]) == 0 ||
        strcmp("pre", argv[i]) == 0 ||
        strcmp("Pre", argv[i]) == 0 ||
        strcmp("pr", argv[i]) == 0 ||
        strcmp("Pr", argv[i]) == 0) {
      command = "preprocess";
      continue;
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
      continue;
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
      continue;
    }
    else if (strcmp("all", argv[i]) == 0 ||
             strcmp("All", argv[i]) == 0) {
      command = "all";
      continue;
    }

    // Check for '--retain' flag, which tells Semblans not to delete outputs from
    // intermediate steps in the pipeline
    else if (strcmp("--retain", argv[i]) == 0 ||
             strcmp("--Retain", argv[i]) == 0 ||
             strcmp("-f", argv[i]) == 0 ||
             strcmp("-F", argv[i]) == 0) {
      retainInterFiles = true;
      continue;
    }
    // Check for '--verbose' flag, which tells Semblans to print a detailed output
    // to the terminal's standard output. This does not affect verbosity of log files.
    // Regardless of terminal output, log files always receive maximum verbosity.
    else if (strcmp("--verbose", argv[i]) == 0 ||
             strcmp("--Verbose", argv[i]) == 0 ||
             strcmp("-v", argv[i]) == 0 ||
             strcmp("-V", argv[i]) == 0) {
      verboseOutput = true;
      continue;
    }
    // If the argument was none of the above, and is the terminating argument,
    // exit, since all subsequent arguments take input
    if (i == argc - 1 && i != nonFlagInd.back()) {
      std::cerr << "\nERROR: Invalid usage. For instructions, call" << std::endl;
      std::cerr << "  semblans --help\n" << std::endl;
      exit(1);
    }
    // Check for config file path
    if (strcmp("--config", argv[i]) == 0 ||
        strcmp("-cfg", argv[i]) == 0) {
      if (strcmp(argv[i + 1] + strlen(argv[i + 1]) - 4, ".ini") == 0 ||
          strcmp(argv[i + 1] + strlen(argv[i + 1]) - 4, ".INI") == 0) {
        pathConfig = argv[i + 1];
        nonFlagInd.push_back(i + 1);
        fs::path pathConfigFile(pathConfig.c_str());
        if (!fs::exists(pathConfigFile)) {
          // ERROR: Config file not found!
          std::cout << "\nERROR: Config file: " << pathConfig << " not found\n" << std::endl;
          exit(1);
        }
      }
      else {
        // ERROR: Config flag invoked, but no config file specified!
        std::cerr << "\nERROR: If using '--config', you must specify config file (.INI)" << std::endl;
        std::cerr << "  (example: --config path/to/config.ini)\n" << std::endl;
        exit(1);
      }
    }
    // Check for '--left' read sequence FASTQ file(s)
    else if (strcmp("--left", argv[i]) == 0 ||
        strcmp("-1", argv[i]) == 0) {
      if (i == argc - 1) {
        std::cerr << "\nERROR: Invalid usage. For instructions, call:\n" << std::endl;
        std::cerr << "  semblans --help\n" << std::endl;
        exit(1);
      }
      leftReads = argv[i + 1];
      nonFlagInd.push_back(i + 1);
    }
    // Check for '--right' read sequence FASTQ file(s)
    else if (strcmp("--right", argv[i]) == 0 ||
             strcmp("-2", argv[i]) == 0) {
      rightReads = argv[i + 1];
      nonFlagInd.push_back(i + 1);
    }
    // Check for '--assembly' transcripts sequence FASTA file
    else if (strcmp("--assembly", argv[i]) == 0 ||
             strcmp("-a", argv[i]) == 0) {
      assembly = argv[i + 1];
      nonFlagInd.push_back(i + 1);
      if (!fs::exists(fs::path(assembly.c_str()))) {
        std::cerr << "ERROR: Assemble '" + assembly + "' does not exist\n" << std::endl;
        exit(1);
      }
    }
    // Check for '--reference-proteome' FASTA file
    else if (strcmp("--reference-proteome", argv[i]) == 0 ||
             strcmp("--ref-proteome", argv[i]) == 0 ||
             strcmp("-rp", argv[i]) == 0) {
      refProt = argv[i + 1];
      nonFlagInd.push_back(i + 1);
      if (!fs::exists(fs::path(refProt.c_str()))) {
        std::cerr << "ERROR: Reference proteome '" + refProt + "' does not exist\n" << std::endl;
        exit(1);
      }
    }
    // Check for '--kraken-db', for specifying which kraken databases to use during foreign
    // sequence filtration
    else if (strcmp("--kraken-db", argv[i]) == 0 ||
             strcmp("-kdb", argv[i]) == 0) {
      kraken2Dbs = argv[i + 1];
      nonFlagInd.push_back(i + 1);
    }
    // Check for '--output-directory', for specifying where outputs should
    // go if no config file is used
    else if (strcmp("--output-directory", argv[i]) == 0 ||
             strcmp("--out-directory", argv[i]) == 0 ||
             strcmp("-o", argv[i]) == 0) {
      outDir = argv[i + 1];
      nonFlagInd.push_back(i + 1);
      if (!fs::exists(outDir.c_str())) {
        std::cerr << "ERROR: Output directory '" + outDir + "' does not exist\n" << std::endl;
        exit(1);
      }
    }
    // Check for '--out', for specifying the name/prefix of the output file
    else if (strcmp("--prefix", argv[i]) == 0 ||
             strcmp("--pre", argv[i]) == 0 ||
             strcmp("-p", argv[i]) == 0) {
      outPrefix = argv[i + 1];
      nonFlagInd.push_back(i + 1);
    }
    // Check for number of threads flag
    else if (strcmp("--threads", argv[i]) == 0 ||
             strcmp("--Threads", argv[i]) == 0 ||
             strcmp("-t", argv[i]) == 0 ||
             strcmp("-T", argv[i]) == 0) {
      if (i != argc - 1) {
        threadStr = argv[i + 1];
        nonFlagInd.push_back(i + 1);
        for (int j = 0; j < strlen(argv[i + 1]); j++) {
          if (!isdigit(threadStr[j])) {
            // ERROR: Bad thread argument
            std::cerr << "ERROR: Invalid argument given for --threads/-t" << std::endl;
            std::cerr << "  (example: --threads 8)\n" << std::endl;
            exit(1);
          }
        }
      }
      else {
        // No thread count given, using 1 thread
        threadStr = "1";
      }
    }

    // Check for ammount of ram to be dedicated
    else if (strcmp("--ram", argv[i]) == 0 ||
             strcmp("--Ram", argv[i]) == 0 ||
             strcmp("-r", argv[i]) == 0 ||
             strcmp("-R", argv[i]) == 0) {

      if (i != argc - 1) {
        ramStr = argv[i + 1];
        nonFlagInd.push_back(i + 1);
        for (int j = 0; j < strlen(argv[i + 1]); j++) {
          if (!isdigit(ramStr[j])) {
            // ERROR: Bad memory argument
            std::cerr << "ERROR: Invalid argument given for --ram/-r" << std::endl;
            std::cerr << "  (example: --ram 10)\n" << std::endl;
            exit(1);
          }
        }
      }
      else {
        // No RAM amount given, using 8 GB
        ramStr = "8";
      }
    }

    // Flag not recognized. Report and exit;
    else {
      for (int k = 0; k < nonFlagInd.size(); k++) {
        if (i == nonFlagInd[k]) {
          argIsFlag = false;
        }
      }
      if (argIsFlag) {
        std::cerr << "\nERROR: Unrecognized flag '" << argv[i] << "'" << std::endl;
        std::cerr << "  For a list of acceptable flags, call" << std::endl;
        std::cerr << "  semblans -h/--help\n" << std::endl;
        exit(1);
      }
    }
  }
}


void checkCommandIsValid(std::string command,
                         std::string leftReads, std::string rightReads,
                         std::string assembly, std::string refProt,
                         std::string outDir, std::string outPrefix,
                         std::string pathConfig, std::string kraken2Dbs) {
  // Check if Semblans command in quick mode (no .INI file) is valid
  if (pathConfig == "null") {
    // Make sure user specified read files (--left/--right) for assembly
    if (leftReads == "null" && rightReads == "null") {
      std::cerr << "\nERROR: User must specify the left/right read file(s) for assembly" << std::endl;
      std::cerr << "  (example: --left reads_1_left.fq,reads_2_left.fq,..." << std::endl;
      std::cerr << "            --right reads_1_right.fq,reads2_right.fq,..." << std::endl;
      exit(1);
    }
    // Make sure user specified an assembly output prefix (--prefix) if they are assembling 
    if (command != "preprocess" && outPrefix == "null") {
      std::cerr << "\nERROR: User must name their assembly with the '--prefix/--pre/-p' flag" << std::endl;
      std::cerr << "  (example: --prefix/--pre/-p assemblyPrefix)" << std::endl;
      exit(1);
    }
    // If user called preprocess, alert user to extraneous command arguments, then proceed
    if (command == "preprocess") {
      if (assembly != "null" || refProt != "null") {
        std::cerr << "\nWARNING: Preprocessing of reads was called, and user specified" << std::endl;
        std::cerr << "an assembly and/or a reference proteome." << std::endl;
        std::cerr << "\n  Neither are necessary for the preprocessing of raw reads," << std::endl;
        std::cerr << "  so they will be ignored" << std::endl;
      }
    }
    // If user called assemble, alert user to extraneous command arguments, then proceed
    if (command == "assemble") {
      if (kraken2Dbs != "null" || refProt != "null") {
        std::cerr << "\nWARNING: Assembly of reads was called, and user specified" << std::endl;
        std::cerr << "a Kraken2 database and/or reference proteome." << std::endl;
        std::cerr << "\n  Neither are necessary for the assembly of reads," << std::endl;
        std::cerr << "  so they will be ignored" << std::endl;
      }
    } 
    // If user called postprocess, make sure they specified an assembly file, reference
    // proteome
    if (command == "postprocess") {
      if (assembly == "null" || refProt == "null") {
        std::cerr << "\nERROR: Postprocessing of transcripts requires an assembly file" << std::endl;
        std::cerr << "and reference proteome for postprocessing" << std::endl;
        std::cerr << "  (example: --assembly/-a transcripts.fa" << std::endl;
        std::cerr << "            --ref-proteome/-rp reference_prot.fa)" << std::endl;
        exit(1);
      }
    }
    // If user is running entire pipeline, make sure they specified a reference proteome
    if (command == "all") {
      if (refProt == "null") {
        std::cerr << "\nERROR: Running the entire pipeline requires a reference proteome" << std::endl;
        std::cerr << "  (example: --ref-proteome/-rp reference_prot.fa)" << std::endl;
        exit(1);
      }
    }
  }
  // Check if Semblans command in advanced mode (with .INI file) is valid
  else {
    if (leftReads != "null" || rightReads != "null" ||
        assembly != "null" || refProt != "null" ||
        outDir != "null" || outPrefix != "null" ||
        kraken2Dbs != "null") {
      std::cerr << "\nERROR: If running in advanced mode (with --config), user should not use the flags:" << std::endl;
      std::cerr << "  --left/-1, --right/-2" << std::endl;
      std::cerr << "  --kraken-db/-kdb" << std::endl;
      std::cerr << "  --assembly/-a" << std::endl;
      std::cerr << "  --ref-proteome/-rp" << std::endl;
      std::cerr << "  --out-directory/-o" << std::endl;
      std::cerr << "These settings are controlled via the .INI configuration file" << std::endl;
      exit(1);
    }
  }
}


// Call the Semblans 'preprocess' binary, which performs all steps up to assembly
int preprocess(std::vector<std::string> sraRuns, bool serialProcess, std::string pathConfig,
               std::string leftReads, std::string rightReads, std::string kraken2Dbs,
               std::string outDir, std::string threadStr, std::string ramStr,
               std::string retain, std::string verbose, std::string entirePipeline,
               std::string logFilePath) {
  
  std::ifstream cfgIniFile;
  std::ofstream cfgIniSub;
  std::string currLine;
  std::string currCfgIniSub;
  std::string currPreCmd;
  size_t pos;
  int result;

  std::string preCmd = SEMBLANS_DIR + "preprocess " + pathConfig + " " +
                       leftReads + " " + rightReads + " " + kraken2Dbs + " " +
                       outDir + " " + threadStr + " " + ramStr +
                       retain + verbose + entirePipeline;

  if (serialProcess) {
    for (auto sraStr : sraRuns) {
      currCfgIniSub = pathConfig.substr(0, pathConfig.rfind(".ini")).insert(pathConfig.rfind("/") + 1, ".") +
                      "." + std::string(fs::path(sraStr.c_str()).stem().c_str()) + ".ini";
      std::replace(currCfgIniSub.begin(), currCfgIniSub.end(), ' ', '_');
      cfgIniFile.open(pathConfig);
      cfgIniSub.open(currCfgIniSub);
      while (std::getline(cfgIniFile, currLine)) {
        for (auto sraStrOther : sraRuns) {
          if (sraStrOther == sraStr) {
            continue;
          }
          pos = currLine.find(sraStrOther);
          if (pos != std::string::npos) {
            currLine.replace(currLine.find(sraStrOther), sraStrOther.length(), "");
          }
        }
        cfgIniSub << currLine << std::endl;
      }
      cfgIniFile.close();
      cfgIniSub.close();
      currPreCmd = SEMBLANS_DIR + "preprocess " + currCfgIniSub + " " +
                   leftReads + " " + rightReads + " " + kraken2Dbs + " " +
                   outDir + " " + threadStr + " " + ramStr + " " + retain +
                   verbose + entirePipeline;

      if (verbose == " true") {
        logOutput("  Running command: " + currPreCmd + "\n", logFilePath);
      }

      result = system(currPreCmd.c_str());
      checkExitSignal(result, logFilePath);
    }
  }
  else {
    if (verbose == " true") {
      logOutput("  Running command: " + preCmd + "\n", logFilePath);
    }
    result = system(preCmd.c_str());
    checkExitSignal(result, logFilePath);
  }
  return result;
}

// Call the Semblans 'assemble' binary, which performs the construction of transcripts from
// short reads
int assemble(std::vector<std::string> sraRuns, bool serialProcess, std::string pathConfig,
             std::string leftReads, std::string rightReads, std::string kraken2Dbs,
             std::string outDir, std::string outPrefix, std::string threadStr, std::string ramStr,
             std::string retain, std::string verbose, std::string entirePipeline,
             std::string logFilePath) {

  std::ifstream cfgIniFile;
  std::ofstream cfgIniSub;
  std::string currLine;
  std::string currCfgIniSub;
  std::string currAssCmd;
  size_t pos;
  int result;

  std::string assCmd = SEMBLANS_DIR + "assemble " + pathConfig + " " +
                       leftReads + " " + rightReads + " " + outDir + " " + outPrefix + " " +
                       threadStr + " " + ramStr + retain + verbose + entirePipeline;

  if (serialProcess) {
    for (auto sraStr : sraRuns) {
      currCfgIniSub = pathConfig.substr(0, pathConfig.rfind(".ini")).insert(pathConfig.rfind("/") + 1, ".") +
                      "." + std::string(fs::path(sraStr.c_str()).stem().c_str()) + ".ini";
      std::replace(currCfgIniSub.begin(), currCfgIniSub.end(), ' ', '_');
      cfgIniFile.open(pathConfig);
      cfgIniSub.open(currCfgIniSub);
      while (std::getline(cfgIniFile, currLine)) {
        for (auto sraStrOther : sraRuns) {
          if (sraStrOther == sraStr) {
            continue;
          }
          pos = currLine.find(sraStrOther);
          if (pos != std::string::npos) {
            currLine.replace(currLine.find(sraStrOther), sraStrOther.length(), "");
          }
        }
        cfgIniSub << currLine << std::endl;
      }
      cfgIniFile.close();
      cfgIniSub.close();

      currAssCmd = SEMBLANS_DIR + "assemble " + currCfgIniSub + " " +
                   leftReads + " " + rightReads + " " + outDir + " " + outPrefix + " " +
                   threadStr + " " + ramStr + retain + verbose + entirePipeline;

      if (verbose == " true") {
        logOutput("  Running command: " + currAssCmd + "\n", logFilePath);
      }
      result = system(currAssCmd.c_str());
      checkExitSignal(result, logFilePath);
    }
  }
  else {
    if (verbose == " true") {
      logOutput("  Running command: " + assCmd + "\n", logFilePath);
    }
    result = system(assCmd.c_str());
    checkExitSignal(result, logFilePath);
  }
  return result;
}

// Call the Semblans 'postprocess' binary, which performs further refining of transcripts
int postprocess(std::vector<std::string> sraRuns, bool serialProcess, std::string pathConfig,
                std::string leftReads, std::string rightReads, std::string assembly,
                std::string refProt, std::string outDir, std::string outPrefix,
                std::string threadStr, std::string ramStr, std::string retain,
                std::string verbose, std::string entirePipeline, std::string logFilePath) {

  std::ifstream cfgIniFile;
  std::ofstream cfgIniSub;
  std::string currLine;
  std::string currCfgIniSub;
  std::string currPostCmd;
  size_t pos;
  int result;

  

  std::string postCmd = SEMBLANS_DIR + "postprocess " + pathConfig + " " +
                        leftReads + " " + rightReads + " " + assembly + " " +
                        refProt + " " + outDir + " " + outPrefix + " " + 
                        threadStr + " " + ramStr + " " +
                        retain + verbose + entirePipeline;

  if (serialProcess) {
    for (auto sraStr : sraRuns) {
      currCfgIniSub = pathConfig.substr(0, pathConfig.rfind(".ini")).insert(pathConfig.rfind("/") + 1, ".") +
                      "." + std::string(fs::path(sraStr.c_str()).stem().c_str()) + ".ini";
      std::replace(currCfgIniSub.begin(), currCfgIniSub.end(), ' ', '_');
      cfgIniFile.open(pathConfig);
      cfgIniSub.open(currCfgIniSub);
      while (std::getline(cfgIniFile, currLine)) {
        for (auto sraStrOther : sraRuns) {
          if (sraStrOther == sraStr) {
            continue;
          }
          pos = currLine.find(sraStrOther);
          if (pos != std::string::npos) {
            currLine.replace(currLine.find(sraStrOther), sraStrOther.length(), "");
          }
        }
        cfgIniSub << currLine << std::endl;
      }
      cfgIniFile.close();
      cfgIniSub.close();

      currPostCmd = SEMBLANS_DIR + "postprocess " + pathConfig + " " +
                    leftReads + " " + rightReads + " " + assembly + " " +
                    refProt + " " + outDir + " " + outPrefix + " " +
                    threadStr + " " + ramStr +
                    retain + verbose + entirePipeline;

      if (verbose == " true") {
        logOutput("  Running command: " + currPostCmd + "\n", logFilePath);
      }
      result = system(currPostCmd.c_str());
      checkExitSignal(result, logFilePath);
    }
  }
  else {
    if (verbose == " true") {
      logOutput("  Running command: " + postCmd + "\n", logFilePath);
    }
    result = system(postCmd.c_str());
    checkExitSignal(result, logFilePath);
  }
  return result;
}


std::string getDirectory(std::string parentDir, std::string dir) {
  std::vector<std::string> matching_files;
  fs::directory_iterator end_itr;
  boost::regex filter(dir);
  for (fs::directory_iterator i(parentDir); i != end_itr; i++) {
    boost::smatch what;
    if (!boost::regex_match(i->path().filename().string(), what, filter)) {
      continue;
    }
    matching_files.push_back(i->path().string());
  }
  return matching_files.front();
}

int main(int argc, char * argv[]) {
  //print_intro();
  INI_MAP cfgIni;
  INI_MAP_ENTRY cfgIniGen;
  std::vector<std::string> sraRuns;
  std::string logFilePath;
  std::string threadStr = "1";
  std::string ramStr = "8";
  int numThreads;
  int ram;
  std::string pathConfig = "null";
  std::string command;
  std::string leftReads = "null";
  std::string rightReads = "null";
  std::string assembly = "null";
  std::string refProt = "null";
  std::string kraken2Dbs = "null";
  std::string outDir = "null";
  std::string outPrefix = "null";
  bool useCfg = true;
  bool serialProcess = false;
  bool retainInterFiles;
  bool verboseOutput;

  std::string preCmd;
  std::string assCmd;
  std::string postCmd;
  if (argc == 1 || 
      (argc == 2 && 
       (strcmp("--help", argv[1]) == 0 ||
        strcmp("-h", argv[1]) == 0))) {
    print_help_base();
    exit(0);
  }
  if (argc == 2 && (strcmp("--version", argv[1]) == 0 ||
                    strcmp("-v", argv[1]) == 0)) {
    std::cout << "Semblans version: v1.0.1" << std::endl;
  }
  else {

    // Parse through flags in commands
    numThreads = 1;
    ram = 8;
    retainInterFiles = false;
    verboseOutput = false;

    // Parse through Semblans argument vector
    parseArgv(argc, argv, command, leftReads, rightReads,
              assembly, refProt, outDir, outPrefix, 
              threadStr, ramStr, pathConfig, kraken2Dbs,
              retainInterFiles, verboseOutput);

    // If no Semblans submodule specified, run all three (preprocess, assemble, postprocess)
    if (command == "") {
      command = "all";
    }

    // Make sure Semblans command entered by user is valid
    checkCommandIsValid(command, leftReads, rightReads, assembly,
                        refProt, outDir, outPrefix, pathConfig, kraken2Dbs);

    // If no config file specified, check for sequence files in Semblans call
    if (pathConfig == "null") {
      useCfg = false;
      if (outDir == "null") {
        system("mkdir semblans_out > /dev/null 2>&1");
        outDir = "semblans_out";
      }
      // If running entire assembly without config file, define name of resulting assembly
      if (command == "all") {
        assembly = outDir + "/assembly/01-Transcript_assembly/" + outPrefix + ".Trinity.fasta";
      }
    }

    // If config file being used, set up workspace for project
    if (useCfg) {
      cfgIni = make_ini_map(pathConfig.c_str());
      cfgIniGen = cfgIni["General"];
      logFilePath = std::string((fs::canonical(fs::path(cfgIniGen["log_file"].c_str()).parent_path()) /
                                 fs::path(cfgIniGen["log_file"].c_str()).filename()).c_str());
    
      fs::create_directory((fs::canonical(fs::path(cfgIniGen["output_directory"].c_str())) /
                            fs::path(cfgIniGen["project_name"].c_str())));
      //std::ofstream logFile(logFilePath, std::ios_base::trunc);

      serialProcess = ini_get_bool(cfgIniGen.at("serial_processing").c_str(), 0);
      if (serialProcess) {
        for (auto sraCode : cfgIni.at("SRA accessions")) {
          sraRuns.push_back(sraCode.first);
        }
        for (auto sraFile : cfgIni.at("Local files")) {
          sraRuns.push_back(sraFile.first);
        }
      }
    }
    // If config file not being used
    else {
      pathConfig = "null";
      
      logFilePath = std::string((fs::canonical(fs::path(outDir.c_str()).parent_path()) /
                                "log.txt").c_str());
    }

    std::ofstream logFile(logFilePath, std::ios_base::trunc);

    // Print Semblans intro in terminal
    print_intro(logFilePath);
 
    std::string verbose;
    std::string retain;
    std::string entirePipeline;
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
    if (command == "all") {
      entirePipeline = " true";
    }
    else {
      entirePipeline = " false";
    }
    int result;



    std::ifstream cfgIniFile;
    std::ofstream cfgIniSub;
    std::string currLine;
    std::string currCfgIniSub;
    std::string currPreCmd;
    std::string currAssCmd;
    std::string currPostCmd;
    size_t pos;

    // Case 1: preprocess
    if (command == "preprocess") {
      logOutput("\nPerforming preprocessing only\n\n", logFilePath);
      preprocess(sraRuns, serialProcess, pathConfig, leftReads, rightReads,
                 kraken2Dbs, outDir, threadStr, ramStr, retain, verbose,
                 entirePipeline, logFilePath);
      exit(0);
    }
    // Case 2: assemble
    if (command == "assemble") {
      logOutput("\nPerforming assembly only\n\n", logFilePath);
      assemble(sraRuns, serialProcess, pathConfig, leftReads, rightReads,
               kraken2Dbs, outDir, outPrefix, threadStr, ramStr, retain, verbose,
               entirePipeline, logFilePath);
      exit(0);
    }
    // Case 3: postprocess
    if (command == "postprocess") {
      logOutput("\nPerforming postprocess only\n\n", logFilePath);
      postprocess(sraRuns, serialProcess, pathConfig, leftReads, rightReads,
                  assembly, refProt, outDir, outPrefix, threadStr, ramStr, retain, verbose,
                  entirePipeline, logFilePath);
      exit(0);
    }
    // Case 4: all three
    if (command == "all") {
      std::string cleanedReadsDir;
      std::vector<std::string> leftReadsCleanedVec;
      std::vector<std::string> rightReadsCleanedVec;
      std::string leftReadsCleaned;
      std::string rightReadsCleaned;
      std::string currLeft;
      std::string currRight;
      std::vector<std::string> assemblies;

      if (pathConfig == "null" && outPrefix == "null") {
        std::cerr << "\nERROR: If not using '--config', user must name their assembly with the '--prefix/--pre' flag" << std::endl;
        std::cerr << "  (example: --prefix/--pre/-p assemblyName" << std::endl;
        exit(1);
      }

      if (serialProcess) {
        for (auto sraStr : sraRuns) {
          currCfgIniSub = pathConfig.substr(0, pathConfig.rfind(".ini")).insert(pathConfig.rfind("/") + 1, ".") +
                          "." + std::string(fs::path(sraStr.c_str()).stem().c_str()) + ".ini";
          
          std::replace(currCfgIniSub.begin(), currCfgIniSub.end(), ' ', '_');
          cfgIniFile.open(pathConfig);
          cfgIniSub.open(currCfgIniSub);
          while (std::getline(cfgIniFile, currLine)) {
            for (auto sraStrOther : sraRuns) {
              if (sraStrOther == sraStr) {
                continue;
              }
              pos = currLine.find(sraStrOther);
              if (pos != std::string::npos) {
                currLine.replace(currLine.find(sraStrOther), sraStrOther.length(), "");
              }
            }
            cfgIniSub << currLine << std::endl;
          }
          cfgIniFile.close();
          cfgIniSub.close();

          std::replace(sraStr.begin(), sraStr.end(), ' ', '\n');

          currPreCmd = SEMBLANS_DIR + "preprocess " + currCfgIniSub + " " +
                       leftReads + " " + rightReads + " " + kraken2Dbs + " " +
                       outDir + " " +
                       std::to_string(numThreads) + " " +
                       std::to_string(ram) + retain + verbose + entirePipeline;

          currAssCmd = SEMBLANS_DIR + "assemble " + currCfgIniSub + " " +
                       leftReads + " " + rightReads + " " + outDir + " " +
                       outPrefix + " " + std::to_string(numThreads) + " " +
                       std::to_string(ram) + retain + verbose + entirePipeline;

          currPostCmd = SEMBLANS_DIR + "postprocess " + currCfgIniSub + " " +
                        leftReads + " " + rightReads + " " + assembly + " " +
                        refProt + " " + outDir + " " + outPrefix + " " +
                        std::to_string(numThreads) + " " +
                        std::to_string(ram) + retain + verbose + entirePipeline;



          logOutput("Performing entire assembly on:\n" + sraStr, logFilePath);
          logOutput("\n ┌───────────────────────────────────────────────────────┐", logFilePath);
          logOutput("\n │         Phase 1: Preprocessing of Short-reads         │", logFilePath);
          logOutput("\n └───────────────────────────────────────────────────────┘\n", logFilePath);
          if (verbose == " true") {
            logOutput("  Running command: " + currPreCmd + "\n", logFilePath);
          }
          result = system(currPreCmd.c_str());
          if (WIFSIGNALED(result) || (result != 0 && WIFEXITED(result) == 1)) {
            checkExitSignal(result, logFilePath);
            system("setterm -cursor on");
            exit(1);
          }
          logOutput("\n ┌───────────────────────────────────────────────────────┐", logFilePath);
          logOutput("\n │        Phase 2: De Novo Assembly of Short-reads       │", logFilePath);
          logOutput("\n └───────────────────────────────────────────────────────┘\n", logFilePath);
          if (verbose == " true") {
            logOutput("  Running command: " + currAssCmd + "\n", logFilePath);
          }
          result = system(currAssCmd.c_str());
          if (WIFSIGNALED(result) || (result != 0 && WIFEXITED(result) == 1)) {
            std::cerr << "\nAssembly exited" << std::endl;
            system("setterm -cursor on");
            exit(1);
          }
          logOutput("\n ┌────────────────────────────────────────────────────────┐", logFilePath);
          logOutput("\n │    Phase 3: Postprocessing of Assembled Transcripts    │", logFilePath);
          logOutput("\n └────────────────────────────────────────────────────────┘\n", logFilePath);
          if (verbose == " true") {
            logOutput("  Running command: " + currPostCmd + "\n", logFilePath);
          }
          result = system(currPostCmd.c_str());
          if (WIFSIGNALED(result) || (result != 0 && WIFEXITED(result) == 1)) {
           std::cerr << "\nPostprocess exited" << std::endl;
           system("setterm -cursor on");
           exit(1);
          }
        }
      }
      else {
        logOutput("Performing entire assembly\n", logFilePath);
        logOutput("\n ┌───────────────────────────────────────────────────────┐", logFilePath);
        logOutput("\n │         Phase 1: Preprocessing of Short-reads         │", logFilePath);
        logOutput("\n └───────────────────────────────────────────────────────┘\n", logFilePath);
        result = preprocess(sraRuns, serialProcess, pathConfig, leftReads, rightReads,
                            kraken2Dbs, outDir, threadStr, ramStr, retain, verbose,
                            entirePipeline, logFilePath);
        if (WIFSIGNALED(result) || (result != 0 && WIFEXITED(result) == 1)) {
          std::cerr << "\nPreprocess exited" << std::endl;
          system("setterm -cursor on");
          exit(1);
        }
        if (pathConfig == "null") {
          cleanedReadsDir = getDirectory(outDir + "/preprocess/", "..-Filter_overrepresented");
          leftReadsCleanedVec = makeCleanedNames(leftReads);
          if (rightReads != "null") {
            rightReadsCleanedVec = makeCleanedNames(rightReads);
          }
          for (int i = 0; i < leftReadsCleanedVec.size(); i++) {
            leftReadsCleaned += cleanedReadsDir + "/" + leftReadsCleanedVec[i];
            if (rightReads != "null") {
              rightReadsCleaned += cleanedReadsDir + "/" + rightReadsCleanedVec[i];
            }
            if (i < leftReadsCleanedVec.size() - 1) {
              leftReadsCleaned += ",";
              rightReadsCleaned += ",";
            }
          }
          if (rightReads == "null") {
            rightReadsCleaned = "null";
          }
        }
        else {
          leftReadsCleaned = "null";
          rightReadsCleaned = "null";
        }

        logOutput("\n ┌───────────────────────────────────────────────────────┐", logFilePath);
        logOutput("\n │        Phase 2: De Novo Assembly of Short-reads       │", logFilePath);
        logOutput("\n └───────────────────────────────────────────────────────┘\n", logFilePath);
        result = assemble(sraRuns, serialProcess, pathConfig, leftReadsCleaned, rightReadsCleaned,
                          kraken2Dbs, outDir, outPrefix, threadStr, ramStr, retain, verbose,
                          entirePipeline, logFilePath);
        if (WIFSIGNALED(result) || (result != 0 && WIFEXITED(result) == 1)) {
          std::cerr << "\nAssembly exited" << std::endl;
          system("setterm -cursor on");
          exit(1);
        }
        logOutput("\n ┌────────────────────────────────────────────────────────┐", logFilePath);
        logOutput("\n │    Phase 3: Postprocessing of Assembled Transcripts    │", logFilePath);
        logOutput("\n └────────────────────────────────────────────────────────┘\n", logFilePath);
        result = postprocess(sraRuns, serialProcess, pathConfig, leftReadsCleaned, rightReadsCleaned,
                             assembly, refProt, outDir, outPrefix, threadStr, ramStr, retain, verbose,
                             entirePipeline, logFilePath);
        if (WIFSIGNALED(result) || (result != 0 && WIFEXITED(result) == 1)) {
          std::cerr << "\nPostprocess exited" << std::endl;
          system("setterm -cursor on");
          exit(1);
        }
        exit(0);
      }
    }
  }
}
