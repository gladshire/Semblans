#include "rem_unfixable.h"


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
    if (currLine1 == "+" || 
       (currLine1.substr(currLine1.length()-5, 5) != "error" &&
        currLine2.substr(currLine2.length()-5, 5) != "error")) {
      outFile1 << currLine1 << std::endl;
      outFile2 << currLine2 << std::endl;
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
    if (currLine == "+" || currLine.substr(currLine.length()-5, 5) != "error") {
      outFile << currLine << std::endl;
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
  int threadNum = stoi(threads);
  std::vector<std::thread> procVec;
  std::vector<std::thread>::iterator start = procVec.begin();
  for (int i = 0; i < sras.size(); i++) {
    /*if (i >= threadNum) {
      std::cout << "Waiting to finish: " << sras[i - threadNum].get_accession() << std::endl;
      procVec.front().join();
      procVec.erase(procVec.begin());
    }*/
    if (sras[i].is_paired()) {
      procVec.push_back(std::thread(rem_unfix_pe, sras[i]));
      std::cout << "Created process for: " << sras[i].get_accession() << std::endl;
    }
    else {
      procVec.push_back(std::thread(rem_unfix_se, sras[i]));
      std::cout << "Created process for: " << sras[i].get_accession() << std::endl;
    }
    if (i == sras.size() - 1) {
      procVec.back().join();
    }
    else {
      procVec.back().detach();
    }
  }
}
