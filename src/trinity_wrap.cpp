#include "trinity_wrap.h"


// Given several SRA run sequence data files, combine into one single file
std::string combine_reads(std::vector<std::string> sraRuns, //std::vector<std::pair<std::string, std::string>> sraRuns,
                          std::string outFileStr,
                          long long int ram_b, std::string logFile) {
  std::string inFileStr1;
  std::string inFileBuffer;
  inFileBuffer.reserve(ram_b);
  std::streamsize s;

  std::ofstream outFile;
  std::ifstream inFile;
  if (fs::exists(fs::path(outFileStr.c_str()))) {
    return outFileStr;
  }
  for (auto sraRun : sraRuns) {
    inFile.open(sraRun.c_str());
    outFile.open(outFileStr, std::ios_base::app);
    while (!inFile.eof()) {
      inFile.read(&inFileBuffer[0], ram_b);
      s = inFile.gcount();
      outFile.write(&inFileBuffer[0], s);
    }
    outFile << std::endl;
    inFile.close();
    //if (sraRun.second != "") {
    //  inFile.open(sraRun.second.c_str());
    //  while (!inFile.eof()) {
    //    inFile.read(&inFileBuffer[0], ram_b);
    //    s = inFile.gcount();
    //    outFile.write(&inFileBuffer[0], s);
    //  }
    //  inFile.close();
    //}
    outFile << std::endl;
    outFile.close();
  }
  return outFileStr;
}

// Perform a de novo assembly on a single SRA's sequence file or file pair with Trinity
void run_trinity(std::pair<std::string, std::string> sraRun, std::string outFile,
                 std::string threads, std::string ram_gb,
                 bool dispOutput, std::string logFile) {
  // Run Trinity for assembly using single SRA
  std::string inFile1 = sraRun.first;
  std::string inFile2 = sraRun.second;
  std::string trin_cmd;
  std::string printOut;
  logOutput("\n    Single assembly chosen", logFile);
  logOutput("\n    Now assembling de-novo transcriptome for:", logFile);
  logOutput("\n      SRA Run:", logFile);
  logOutput("\n      " + std::string(fs::path(sraRun.first).stem().c_str()), logFile);
  if (sraRun.second != "") {
    logOutput("\n      " + std::string(fs::path(sraRun.second).stem().c_str()), logFile);
  }
  logOutput("\n", logFile);
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  if (inFile2 != "") {
    trin_cmd = PATH_TRINITY + " --seqType fq" +
               " --left " + inFile1 + " --right " + inFile2 +
               " --max_memory " + ram_gb + "G " + "--CPU " + threads +
               " --bflyCalculateCPU" + " --full_cleanup" +
               " --no_normalize_reads" + " --output " + outFile + printOut;
  }
  else {
    trin_cmd = PATH_TRINITY + " --seqType fq" +
               " --single " + inFile1 +
               " --max_memory " + ram_gb + "G " + "--CPU " + threads +
               " --bflyCalculateCPU" + " --full_cleanup" +
               " --no_normalize_reads" + " --output " + outFile + printOut;
  }

  if (!dispOutput) {
    procRunning = true;
    std::thread assembleThread(progressAnim, "    ", logFile);
    result = system(trin_cmd.c_str());
    procRunning = false;
    assembleThread.join();
  }
  else {
    result = system(trin_cmd.c_str());
  }
  if (WIFSIGNALED(result)) {
    system("setterm -cursor on");
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
  std::rename((outFile + ".Trinity.fasta").c_str(), outFile.c_str());
  std::rename((outFile + ".Trinity.fasta.gene_trans_map").c_str(),
              (outFile + ".gene_trans_map").c_str());
  replaceChar(logFile, '\r', '\n');
}

// Perform a de novo assembly using multiple SRAs' sequence data with Trinity
void run_trinity_comb(std::vector<std::pair<std::string, std::string>> sraRuns,
                      std::string outFile,
                      std::string threads, std::string ram_gb,
                      bool dispOutput, std::string logFile) {
  
  // Run Trinity for assembly using multiple SRA runs
  logOutput("\n  Combined assembly chosen\n\n  Preparing reads for assembly\n", logFile);
  long long int ram_b = (long long int)stoi(ram_gb) * 1000000000;

  // Prepare combined SRA file for multiple assembly
  // Determine whether data is mixture of single and paired end
  bool isMixed = false;
  bool combPaired;
  if (sraRuns[0].second != "") combPaired = true;
  for (auto sraRun : sraRuns) {
    if (combPaired && sraRun.second == "") {
      isMixed = true;
    }
    if (!combPaired && sraRun.second != "") {
      isMixed = true;
    }
  }
  std::string inFile1;
  std::string inFile2;
  std::string outComb1;
  std::string outComb2;
  std::string outFileComb1;
  std::string outFileComb2;

  if (!isMixed && combPaired) {
    outComb1 = std::string(fs::path(outFile.c_str()).stem().c_str()) + "_1.comb.fastq";
    outComb2 = std::string(fs::path(outFile.c_str()).stem().c_str()) + "_2.comb.fastq";
    outFileComb1 = std::string(fs::path(outFile.c_str()).parent_path().c_str()) + "/" + outComb1;
    outFileComb2 = std::string(fs::path(outFile.c_str()).parent_path().c_str()) + "/" + outComb2;

    std::vector<std::string> srasLeft;
    std::vector<std::string> srasRight;
    for (auto sraRun : sraRuns) {
      srasLeft.push_back(sraRun.first);
      srasRight.push_back(sraRun.second);
    }

    if (!dispOutput) {
      procRunning = true;
      std::thread combProgThread(progressAnim, "    ", logFile);
      inFile1 = combine_reads(srasLeft, outFileComb1, ram_b, logFile);
      inFile2 = combine_reads(srasRight, outFileComb2, ram_b, logFile);
      procRunning = false;
      combProgThread.join();
    }
    else {
      inFile1 = combine_reads(srasLeft, outFileComb1, ram_b, logFile);
      inFile2 = combine_reads(srasRight, outFileComb2, ram_b, logFile);
    }
  }
  else if (!isMixed && !combPaired) {
    outComb1 = std::string(fs::path(outFile.c_str()).stem().c_str()) + ".comb.fastq";
    outFileComb1 = std::string(fs::path(outFile.c_str()).parent_path().c_str()) + "/" +
                               outComb1;

    std::vector<std::string> srasSingle;
    for (auto sraRun : sraRuns) {
      srasSingle.push_back(sraRun.first);
    }

    if (!dispOutput) {
      procRunning = true;
      std::thread combProgThread(progressAnim, "    ", logFile);
      inFile1 = combine_reads(srasSingle, outFileComb1, ram_b, logFile);
      procRunning = false;
      combProgThread.join();
    }
    else {
      inFile1 = combine_reads(srasSingle, outFileComb1, ram_b, logFile);
    }
  }

  // Summarize Trinity assembly job 
  logOutput("    Now assembling de-novo transcriptome for:\n", logFile);
  for (auto sraRun : sraRuns) {
    logOutput("\n      SRA Run:", logFile);
    logOutput("\n      " + std::string(fs::path(sraRun.first).stem().c_str()), logFile);
    if (sraRun.second != "") {
      logOutput("\n      " + std::string(fs::path(sraRun.second).stem().c_str()), logFile);
    }
    logOutput("\n", logFile);
  }
  
  std::string trin_cmd;
  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }

  int result;
  if (!isMixed) {
    if (combPaired) {
      run_trinity(std::pair<std::string, std::string>(outFileComb1, outFileComb2),
                  outFile, threads, ram_gb, dispOutput, logFile);
    }
    else {
      run_trinity(std::pair<std::string, std::string>(outFileComb1, ""),
                  outFile, threads, ram_gb, dispOutput, logFile);
    }
  }
  else {
    logOutput("    ERROR: Semblans does not support mixing paired and unpaired reads", logFile);
    system("setterm -cursor on");
    exit(1);
    /*
    trin_cmd = PATH_TRINITY + " --seqType fq" + " --single " + inFile1 + " --max_memory " +
               ram_gb + "G " + "--CPU " + threads + " --bflyCalculateCPU" + " --full_cleanup" +
               " --no_normalize_reads" + " --run_as_paired" + " --output " + outFile + printOut; 
  
    // Run Trinity assembly
    procRunning = true;
    std::thread trinProgThread(progressAnim, "    ", logFile);
    result = system(trin_cmd.c_str());
    procRunning = false;
    trinProgThread.join();

    if (WIFSIGNALED(result)) {
      system("setterm -cursor on");
      logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
      exit(1);
    }
    std::rename((outFile + ".Trinity.fasta").c_str(), outFile.c_str());
    std::rename((outFile + ".Trinity.fasta.gene_trans_map").c_str(),
              (outFile + ".gene_trans_map").c_str());
    replaceChar(logFile, '\r', '\n');
    */
  }
}
