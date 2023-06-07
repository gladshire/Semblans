#include "sra_toolkit.h"


std::vector<SRA> get_sras(const INI_MAP &iniFile, bool compressFiles) {
  std::vector<SRA> sras;
  if (!iniFile.at("SRA accessions").empty()) {
    for (auto sra : iniFile.at("SRA accessions")) {
      sras.push_back(SRA(sra.first, iniFile, compressFiles));
    }
  }
  return sras;
}

void prefetch_sra(SRA sra, bool dispOutput, std::string logFile) {
  std::string outDir(sra.get_sra_path_raw().first.parent_path().c_str());
  std::string prefetchFlag = " --max-size u -O ";
  std::string sraAccession = sra.get_accession();
  int result;
  std::string prefetchCmd = PATH_PREFETCH + " " + sraAccession + prefetchFlag + outDir;
  if (dispOutput) {
    prefetchCmd += (" 2>&1 | tee -a " + logFile);
  }
  else {
    prefetchCmd += (" >>" + logFile + " 2>&1");
  }
  result = system(prefetchCmd.c_str());
  if (WIFSIGNALED(result)) {
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}


void fasterq_sra(SRA sra, std::string threads, bool dispOutput,
                 bool compressOutput, std::string logFile) {
  std::string prefetchDir(sra.get_sra_path_raw().first.parent_path().c_str());
  std::string outFile;
  std::string fasterqFlag = " -e " + threads + " -t " + prefetchDir + " " + prefetchDir;
  std::string sraAccession = sra.get_accession();
  int result;
  outFile = sra.make_file_str();
  fs::path currDir = fs::current_path();
  fs::current_path(fs::path(prefetchDir.c_str()));
  std::string fasterqCmd;
  if (compressOutput) {
    fasterqCmd = "( " + PATH_FASTERQ + " " + fasterqFlag + "/" + sraAccession +
                 " --split-spot -Z | awk \'{" +
                 "if ((NR-1) % 8 < 4) {print | \"" + PATH_PIGZ + " --fast -p " + threads + " > " + outFile + "_1.fastq.gz\"} " +
                 "else {print | \"" + PATH_PIGZ + " --fast -p " + threads + " > " + outFile + "_2.fastq.gz\"} }\' )";
  }
  else {
    fasterqCmd = PATH_FASTERQ + " " + fasterqFlag + "/" + sraAccession +
                 " -o " + outFile;
  }
  if (dispOutput) {
    fasterqCmd += (" 2>&1 | tee -a " + logFile);
  }
  else {
    fasterqCmd += (" >>" + logFile + " 2>&1");
  }

  result = system(fasterqCmd.c_str());
  fs::current_path(currDir);
  if (WIFSIGNALED(result)) {
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}


void align_file_buffer(std::ifstream & inFile1, std::ifstream & inFile2,
                       char * inFile1Data, char * inFile2Data,
                       std::streamsize & s1, std::streamsize & s2) {
  std::string read;
  std::string readMatch;
  if (!inFile1.eof() && !inFile2.eof()) {
    while ((inFile1.peek() != '@' && inFile1.peek() != '>') &&
           (inFile2.peek() != '@' && inFile2.peek() != '>')) {
      if (inFile1.peek() != '@' && inFile1.peek() != '>') {
        inFile1.unget();
        inFile1Data[s1 - 1] = '\0';
        s1--;
      }
      if (inFile2.peek() != '@' && inFile2.peek() != '>') {
        inFile2.unget();
        inFile2Data[s2 - 1] = '\0';
        s2--;
      }
    }
    if (inFile1.peek() == '@' || inFile1.peek() == '>') {
      inFile1.get();
      inFile1 >> readMatch;
      while (inFile1.peek() != '@' && inFile1.peek() != '>') {
        inFile1.unget();
      }
      while (inFile2.peek() != '@' && inFile2.peek() != '>') {
        inFile2.unget();
        inFile2Data[s2 - 1] = '\0';
        s2--;
      }
      inFile2.get();
      inFile2 >> read;
      while (inFile2.peek() != '@' && inFile2.peek() != '>') {
        inFile2.unget();
      }
      if (read > readMatch) {
        inFile2.unget();
        while (read.compare(readMatch) != 0) {
          while (inFile2.peek() != '@' && inFile2.peek() != '>') {
            inFile2.unget();
            inFile2Data[s2 - 1] = '\0';
            s2--;
          }
          inFile2.get();
          inFile2 >> read;
          while (inFile2.peek() != '@' && inFile2.peek() != '>') {
            inFile2.unget();
          }
          inFile2.unget();
          inFile2Data[s2 - 1] = '\0';
          s2--;
        }
        inFile1.get();
      }
      else if (read < readMatch) {
        inFile1.unget();
        while (read.compare(readMatch) != 0) {
          while (inFile1.peek() != '@' && inFile1.peek() != '>') {
            inFile1.unget();
            inFile1Data[s1 - 1] = '\0';
            s1--;
          }
          inFile1.get();
          inFile1 >> readMatch;
          while (inFile1.peek() != '@' && inFile1.peek() != '>') {
            inFile1.unget();
          }
          inFile1.unget();
          inFile1Data[s1 - 1] = '\0';
          s1--;
        }
        inFile1.get();
      }
      else {
        // Buffer in position -- do nothing
      }
    }
    else {
      inFile2.get();
      inFile2 >> readMatch;
      while (inFile2.peek() != '@' && inFile2.peek() != '>') {
        inFile2.unget();
      }
      while (inFile1.peek() != '@' && inFile1.peek() != '>') {
        inFile1.unget();
        inFile1Data[s1 - 1] = '\0';
        s1--;
      }
      inFile1.get();
      inFile1 >> read;
      while (inFile1.peek() != '@' && inFile1.peek() != '>') {
        inFile1.unget();
      }
      if (read > readMatch) {
        inFile1.unget();
        while (read.compare(readMatch) != 0) {
          while (inFile1.peek() != '@' && inFile1.peek() != '>') {
            inFile1.unget();
            inFile1Data[s1 - 1] = '\0';
            s1--;
          }
          inFile1.get();
          inFile1 >> read;
          while (inFile1.peek() != '@' && inFile1.peek() != '>') {
            inFile1.unget();
          }
          inFile1.unget();
          inFile1Data[s1 - 1] = '\0';
          s1--;
        }
        inFile1.get();
      }
      else if (read < readMatch) {
        inFile2.unget();
        while (read.compare(readMatch) != 0) {
          while (inFile2.peek() != '@' && inFile2.peek() != '>') {
            inFile2.unget();
            inFile2Data[s2 - 1] = '\0';
            s2--;
          }
          inFile2.get();
          inFile2 >> readMatch;
          while (inFile2.peek() != '@' && inFile2.peek() != '>') {
            inFile2.unget();
          }
          inFile2.unget();
          inFile2Data[s2 - 1] = '\0';
          s2--;
        }
        inFile2.get();
      }
      else {
        // Buffer in position -- do nothing
      }
    }
  }
}
