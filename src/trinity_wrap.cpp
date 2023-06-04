#include "trinity_wrap.h"

// TODO: Change condition by which multi-SRA assembly is determined
//   NOW: SRAs with same organism name are combined
//     - Problem if local files, without organism names are used

std::vector<SRA> get_sra_to_combine(std::vector<SRA> sras, std::string org_name) {
  std::vector<SRA> sra_comb;
  for (auto &sra : sras) {
    // TODO: FIX THIS; SHOULD NOT JUST BE BY ORGANISM NAME
    if (sra.get_org_name() == org_name) {
      sra_comb.push_back(sra);
    }
  }
  return sra_comb;
}

std::string combine_reads(std::vector<std::pair<std::string, std::string>> sraRuns,
                          std::string outFileStr,
                          long long int ram_b, std::string logFile) {
  std::string inFileStr1;
  std::string inFileBuffer;
  inFileBuffer.reserve(ram_b);
  std::streamsize s;

  std::ofstream outFile;
  std::ifstream inFile;
  logOutput("Combined assembly chosen\nPreparing reads for assembly ...", logFile);
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

void run_trinity(std::pair<std::string, std::string> sraRun, std::string outFile,
                 std::string threads, std::string ram_gb,
                 bool dispOutput, std::string logFile) {
  // Run Trinity for assembly using single SRA

  // Summarize Trinity assembly job
  logOutput("Now assembling de-novo transcriptome for:", logFile);
  logOutput("  " + sraRun.first, logFile);
  if (sraRun.second != "") {
    logOutput("  " + sraRun.second, logFile);
  }
  std::string inFile1 = sraRun.first;
  std::string inFile2 = sraRun.second;
  std::string trin_cmd;
  std::string printOut;
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
  result = system(trin_cmd.c_str());
  if (WIFSIGNALED(result)) {
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
  std::rename((outFile + ".Trinity.fasta").c_str(), outFile.c_str());
}

void run_trinity_comb(std::vector<std::pair<std::string, std::string>> sraRuns,
                      std::string outFile,
                      std::string threads, std::string ram_gb,
                      bool dispOutput, std::string logFile) {
  // Run Trinity for assembly using multiple SRAs

  // Summarize Trinity assembly job
  logOutput("Now assembling de-novo transcriptome for:", logFile);
  for (auto sraRun : sraRuns) {
    logOutput("\n  SRA Run:", logFile);
    logOutput("\n  " + sraRun.first, logFile);
    if (sraRun.second != "") {
      logOutput("\n  " + sraRun.second, logFile);
    }
    logOutput("\n", logFile);
  }
  long long int ram_b = (long long int)stoi(ram_gb) * 1000000000;
  std::string outComb(std::string(fs::path(outFile.c_str()).stem().c_str()) + ".comb.fastq");
  std::string outFileComb(std::string(fs::path(outFile.c_str()).parent_path().c_str()) + "/" +
                          outComb);

  std::string inFile = combine_reads(sraRuns, outFileComb, ram_b, logFile);
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
  result = system(trin_cmd.c_str());

  if (WIFSIGNALED(result)) {
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
  std::rename((outFile + ".Trinity.fasta").c_str(), outFile.c_str());
}
