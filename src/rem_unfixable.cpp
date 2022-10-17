#include "rem_unfixable.h"

std::vector<std::thread> procVec;
std::mutex threadMutex;

void rem_unfix_pe(SRA sra) {
  std::string inFile1Str(sra.get_sra_path_corr().first.c_str());
  std::string inFile2Str(sra.get_sra_path_corr().second.c_str());
  std::string outFile1Str(std::string(sra.get_sra_path_corr().first.replace_extension("fix.fq").c_str()));
  std::string outFile2Str(std::string(sra.get_sra_path_corr().second.replace_extension("fix.fq").c_str()));
  std::ifstream inFile1(inFile1Str);
  std::ifstream inFile2(inFile2Str);
  std::ofstream outFile1(outFile1Str);
  std::ofstream outFile2(outFile2Str);

  std::string currLine1;
  std::string currLine2;

  if (!inFile1.is_open() || !inFile2.is_open()) {
    std::cout << "Cannot open file(s)" << std::endl;
    return;
  }
  while (getline(inFile1, currLine1) && getline(inFile2, currLine2)) {
  /*char buffer1[1024];
  char buffer2[1024];
  char * buffer1L;
  char * buffer2L;
  char * nl1c;
  char * nl2c;
  char * nl1l;
  char * nl2l;

  while (!inFile1.eof() && !inFile2.eof()) {
    inFile1.read(buffer1, buffer1.sizeof());
    inFile2.read(buffer2, buffer2.sizeof());
    std::streamsize s1 = inFile1.gcount();
    std::streamsize s2 = inFile2.gcount();
    nl1c = strtok(buffer1, '\n');
    nl2c = strtok(buffer2, '\n');
    nl1l = buffer1;
    nl2l = buffer2;
    while (nlc1 != NULL || nlc2 != NULL) {
      if (strcmp(nl1c - 5, "error") == 0 ||
          strcmp(nl2c - 5, "error") == 0) {
        if (nl1c - buffer1 < nl2c - buffer2) {
          
        }
      }
   */   
    if (currLine1.size() == 1 || 
       (currLine1.compare(currLine1.length()-5, 5, "error") != 0 &&
        currLine2.compare(currLine2.length()-5, 5, "error") != 0)) {
      outFile1 << currLine1 << '\n';
      outFile2 << currLine2 << '\n';
    }
    else {
      for (int i = 0; i < 3; i++) {
        getline(inFile1, currLine1);
        getline(inFile2, currLine2);
      }
    }
  }
  inFile1.close();
  inFile2.close();
  outFile1.close();
  outFile2.close();
}


void rem_unfix_se(SRA sra) {
  std::string inFileStr(sra.get_sra_path_corr().first.c_str());
  std::string outFileStr(std::string(sra.get_sra_path_corr().first.replace_extension("fix.fq").c_str()));

  std::ifstream inFile(inFileStr);
  std::ofstream outFile(outFileStr);

  std::string currLine;

  if (!inFile.is_open()) {
    std::cout << "Cannot open file(s)" << std::endl;
  }
  while (getline(inFile, currLine)) {
    if (currLine.size() == 1 ||
        currLine.compare(currLine.length()-5, 5, "error") != 0) {
      outFile << currLine << '\n';
    }
    else {
      for (int i = 0; i < 3; i++) {
        getline(inFile, currLine);
      }
    }
  }
  inFile.close();
  outFile.close();
}


void rem_unfix_bulk(std::vector<SRA> sras, std::string threads) {
  std::cout << "\nRemoving unfixable reads for:\n" << std::endl;
  summarize_all_sras(sras);

  int threadNum = stoi(threads);
  threadPool fileCorrPool;
  fileCorrPool.start(threadNum);
  for (auto sra : sras) {
    if (fs::exists(sra.get_sra_path_corr_fix().first.c_str())) {
      std::cout << "Fixed version found for: " << sra.get_accession() << std::endl;
      continue;
    }
    if (sra.is_paired()) {
      fileCorrPool.queueJob([sra] {rem_unfix_pe(sra);});
    }
    else {
      fileCorrPool.queueJob([sra] {rem_unfix_se(sra);});
    }
  }
  fileCorrPool.stop();
}
