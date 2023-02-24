#include "sra_toolkit.h"


std::vector<SRA> get_sras(const INI_MAP &iniFile) {
  std::vector<SRA> sras;
  if (!iniFile.at("SRA accessions").empty()) {
    for (auto sra : iniFile.at("SRA accessions")) {
      sras.push_back(SRA(sra.first, iniFile));
    }
  }
  return sras;
}

void prefetch_sra(std::vector<SRA> sras) {
  std::string outDir(sras[0].get_sra_path_raw().first.parent_path().native().c_str());
  std::string sra_accessions = "";
  std::string prefetchFlag = " --max-size u -p -O ";
  int result;
  for (auto sra : sras) {
    if (fs::exists(fs::path(std::string(sra.get_sra_path_raw().first.parent_path().c_str()) +
                   "/" + sra.get_accession()))) {
      std::cout << "Prefetch found for: " << sra.get_accession() << std::endl;
      continue;
    }
    sra_accessions += (sra.get_accession() + " ");
  }
  if (sra_accessions != "") {
    result = system((PATH_PREFETCH + " " + sra_accessions + prefetchFlag + outDir).c_str());
    if (WIFSIGNALED(result)) {
      std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
      exit(1);
    }
  }
}

void fasterq_sra(std::vector<SRA> sras, std::string threads) {
  std::string prefetchDir(sras[0].get_sra_path_raw().first.parent_path().native().c_str());
  std::string outFile;
  std::string fasterqFlag = " -p -e " + threads + " -t " + prefetchDir + " " + prefetchDir;
  int result;
  for (auto sra : sras) {
    if (fs::exists(sra.get_sra_path_raw().first)) {
      std::cout << "Raw reads found for: " << sra.get_accession() << std::endl;
      continue;
    }
    outFile = sra.make_file_str();
    fs::path currDir = fs::current_path();
    fs::current_path(fs::path(prefetchDir.c_str()));
    result = system((PATH_FASTERQ + " " + fasterqFlag + "/" + sra.get_accession() +
                     " -o " + outFile).c_str());
    fs::current_path(currDir);
    if (WIFSIGNALED(result)) {
      std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
      exit(1);
    }
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
