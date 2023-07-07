#include "sra_toolkit.h"


// Given a Semblans config INI file, construct and extract SRA objects
// based on accession numbers
std::vector<SRA> get_sras(const INI_MAP &iniFile, bool compressFiles) {
  std::vector<SRA> sras;
  if (!iniFile.at("SRA accessions").empty()) {
    for (auto sra : iniFile.at("SRA accessions")) {
      sras.push_back(SRA(sra.first, iniFile, compressFiles));
    }
  }
  return sras;
}

// Given an SRA run object, prefetch the raw pre-sequence data using SRA-toolkit
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
    system("setterm -cursor on");
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}

// Given an SRA run object, dump prefetched data to a FASTQ file 
void fasterq_sra(SRA sra, std::string threads, bool dispOutput,
                 bool compressOutput, std::string logFile) {
  std::string prefetchDir(sra.get_sra_path_raw().first.parent_path().c_str());
  std::string outFile;
  std::string fasterqFlag = " -e " + threads + " -t ./";
  std::string sraAccession = sra.get_accession();
  int result;
  outFile = sra.make_file_str();
  if (!sra.is_paired()) {
    outFile += ".fastq";
  }
  fs::path currDir = fs::current_path();
  fs::current_path(fs::path(prefetchDir.c_str()));
  std::string fasterqCmd;
  if (compressOutput) {
    fasterqCmd = "( " + PATH_FASTERQ + " " + sraAccession + fasterqFlag +
                 " --split-spot -Z | awk \'{" +
                 "if ((NR-1) % 8 < 4) {print | \"" + PATH_PIGZ + " --fast -p " + threads + " > " + outFile + "_1.fastq.gz\"} " +
                 "else {print | \"" + PATH_PIGZ + " --fast -p " + threads + " > " + outFile + "_2.fastq.gz\"} }\' )";
  }
  else {
    fasterqCmd = "(" + PATH_FASTERQ + " " + sraAccession + fasterqFlag +
                 " -o ./" + outFile + " )";
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
    system("setterm -cursor on");
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}

// Align in file streams of paired-end sequence data to ensure both data buffers
// contain the same reads
void align_file_buffer(std::istream & inStream1, std::istream & inStream2,
                       char * inStream1Data, char * inStream2Data,
                       std::streamsize & s1, std::streamsize & s2) {
  std::string read;
  std::string readMatch;
  if (!inStream1.eof() && !inStream2.eof()) {
    while ((inStream1.peek() != '@' && inStream1.peek() != '>') &&
           (inStream2.peek() != '@' && inStream2.peek() != '>')) {
      if (inStream1.peek() != '@' && inStream1.peek() != '>') {
        inStream1.unget();
        inStream1Data[s1 - 1] = '\0';
        s1--;
      }
      if (inStream2.peek() != '@' && inStream2.peek() != '>') {
        inStream2.unget();
        inStream2Data[s2 - 1] = '\0';
        s2--;
      }
    }
    if (inStream1.peek() == '@' || inStream1.peek() == '>') {
      inStream1.get();
      inStream1 >> readMatch;
      while (inStream1.peek() != '@' && inStream1.peek() != '>') {
        inStream1.unget();
      }
      while (inStream2.peek() != '@' && inStream2.peek() != '>') {
        inStream2.unget();
        inStream2Data[s2 - 1] = '\0';
        s2--;
      }
      inStream2.get();
      inStream2 >> read;
      while (inStream2.peek() != '@' && inStream2.peek() != '>') {
        inStream2.unget();
      }
      if (read > readMatch) {
        inStream2.unget();
        while (read.compare(readMatch) != 0) {
          while (inStream2.peek() != '@' && inStream2.peek() != '>') {
            inStream2.unget();
            inStream2Data[s2 - 1] = '\0';
            s2--;
          }
          inStream2.get();
          inStream2 >> read;
          while (inStream2.peek() != '@' && inStream2.peek() != '>') {
            inStream2.unget();
          }
          inStream2.unget();
          inStream2Data[s2 - 1] = '\0';
          s2--;
        }
        inStream1.get();
      }
      else if (read < readMatch) {
        inStream1.unget();
        while (read.compare(readMatch) != 0) {
          while (inStream1.peek() != '@' && inStream1.peek() != '>') {
            inStream1.unget();
            inStream1Data[s1 - 1] = '\0';
            s1--;
          }
          inStream1.get();
          inStream1 >> readMatch;
          while (inStream1.peek() != '@' && inStream1.peek() != '>') {
            inStream1.unget();
          }
          inStream1.unget();
          inStream1Data[s1 - 1] = '\0';
          s1--;
        }
        inStream1.get();
      }
      else {
        // Buffer in position -- do nothing
      }
    }
    else {
      inStream2.get();
      inStream2 >> readMatch;
      while (inStream2.peek() != '@' && inStream2.peek() != '>') {
        inStream2.unget();
      }
      while (inStream1.peek() != '@' && inStream1.peek() != '>') {
        inStream1.unget();
        inStream1Data[s1 - 1] = '\0';
        s1--;
      }
      inStream1.get();
      inStream1 >> read;
      while (inStream1.peek() != '@' && inStream1.peek() != '>') {
        inStream1.unget();
      }
      if (read > readMatch) {
        inStream1.unget();
        while (read.compare(readMatch) != 0) {
          while (inStream1.peek() != '@' && inStream1.peek() != '>') {
            inStream1.unget();
            inStream1Data[s1 - 1] = '\0';
            s1--;
          }
          inStream1.get();
          inStream1 >> read;
          while (inStream1.peek() != '@' && inStream1.peek() != '>') {
            inStream1.unget();
          }
          inStream1.unget();
          inStream1Data[s1 - 1] = '\0';
          s1--;
        }
        inStream1.get();
      }
      else if (read < readMatch) {
        inStream2.unget();
        while (read.compare(readMatch) != 0) {
          while (inStream2.peek() != '@' && inStream2.peek() != '>') {
            inStream2.unget();
            inStream2Data[s2 - 1] = '\0';
            s2--;
          }
          inStream2.get();
          inStream2 >> readMatch;
          while (inStream2.peek() != '@' && inStream2.peek() != '>') {
            inStream2.unget();
          }
          inStream2.unget();
          inStream2Data[s2 - 1] = '\0';
          s2--;
        }
        inStream2.get();
      }
      else {
        // Buffer in position -- do nothing
      }
    }
  }
}
