#include "trinity_wrap.h"


// Given several SRA run sequence data files, combine into one single file
std::string combine_reads(std::vector<std::pair<std::string, std::string>> sraRuns,
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
    inFile.open(sraRun.first.c_str());
    outFile.open(outFileStr, std::ios_base::app);
    while (!inFile.eof()) {
      inFile.read(&inFileBuffer[0], ram_b);
      s = inFile.gcount();
      outFile.write(&inFileBuffer[0], s);
    }
    outFile << std::endl;
    inFile.close();
    if (sraRun.second != "") {
      inFile.open(sraRun.second.c_str());
      while (!inFile.eof()) {
        inFile.read(&inFileBuffer[0], ram_b);
        s = inFile.gcount();
        outFile.write(&inFileBuffer[0], s);
      }
      inFile.close();
    }
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
  logOutput("    Single assembly chosen", logFile);
  logOutput("\n    Now assembling de-novo transcriptome for:\n", logFile);
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
}

// Perform a de novo assembly using multiple SRAs' sequence data with Trinity
void run_trinity_comb(std::vector<std::pair<std::string, std::string>> sraRuns,
                      std::string outFile,
                      std::string threads, std::string ram_gb,
                      bool dispOutput, std::string logFile) {
  
  // Run Trinity for assembly using multiple SRA runs
  logOutput("\n  Combined assembly chosen\n\n  Preparing reads for assembly\n", logFile);
  long long int ram_b = (long long int)stoi(ram_gb) * 1000000000;
  std::string outComb(std::string(fs::path(outFile.c_str()).stem().c_str()) + ".comb.fastq");
  std::string outFileComb(std::string(fs::path(outFile.c_str()).parent_path().c_str()) + "/" +
                          outComb);
  // Prepare combined SRA file for multiple assembly
  std::string inFile;
  if (!dispOutput) {
    procRunning = true;
    std::thread combProgThread(progressAnim, "    ", logFile);
    inFile = combine_reads(sraRuns, outFileComb, ram_b, logFile);
    procRunning = false;
    combProgThread.join();
  }
  else {
    inFile = combine_reads(sraRuns, outFileComb, ram_b, logFile);
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
  trin_cmd = PATH_TRINITY + " --seqType fq" + " --single " + inFile + " --max_memory " +
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
}
