#include "rem_unfixable.h"

std::vector<std::thread> procVec;
std::mutex threadMutex;

void rem_unfix_pe(SRA sra, long long int ram_b) {
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

  /*long long int ram_b_per_file = ram_b / 2;

  inFile1.seekg(0, inFile1.end);
  long int lenFile1 = inFile1.tellg();
  inFile1.seekg(0, inFile1.beg);

  inFile2.seekg(0, inFile2.end);
  long int lenFile2 = inFile2.tellg();
  inFile2.seekg(0, inFile2.beg);

  std::string inFile1Data;
  std::string inFile2Data;

  std::streamsize s1;
  std::streamsize s2;

  inFile1Data.reserve(ram_b_per_file);
  inFile2Data.reserve(ram_b_per_file);

  size_t currPos1 = 0;
  size_t currPos2 = 0;
  size_t nlPos1;
  size_t nlPos2;

  while(!inFile1.eof() && !inFile2.eof()) {
    fread(inFile1Data, ram_b_per_file, inFile1);
    fread(inFile2Data, ram_b_per_file, inFile2);
    s1 = inFile1.gcount();
    s2 = inFile2.gcount();

    // Unstream end of buffer to just before last '@' or '>' character in inFile1, inFile2
    // Remove characters from inFile1Data, inFile2Data appropriately
    //
    // Purpose: removes edge cases, allows for 4-line processing, increasing
    //          parse time by: 4RT() - T(unget)

  }*/

  while (getline(inFile1, currLine1) && getline(inFile2, currLine2)) {
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


void rem_unfix_se(SRA sra, long long int ram_b) {
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


void rem_unfix_bulk(std::vector<SRA> sras, std::string threads, std::string ram_gb) {
  std::cout << "\nRemoving unfixable reads for:\n" << std::endl;
  summarize_all_sras(sras);

  int threadNum = stoi(threads);
  long long int ram_b = stoi(ram_gb) * 1000000000 - 1000;
  long long int ram_b_per_thread = ram_b / threadNum;

  threadPool fileCorrPool;
  fileCorrPool.start(threadNum);
  for (auto sra : sras) {
    if (fs::exists(sra.get_sra_path_corr_fix().first.c_str())) {
      std::cout << "Fixed version found for: " << sra.get_accession() << std::endl;
      continue;
    }
    if (sra.is_paired()) {
      fileCorrPool.queueJob([sra, ram_b_per_thread] {rem_unfix_pe(sra, ram_b_per_thread);});
    }
    else {
      fileCorrPool.queueJob([sra, ram_b_per_thread] {rem_unfix_se(sra, ram_b_per_thread);});
    }
  }
  fileCorrPool.stop();
}
