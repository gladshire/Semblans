#include "rem_unfixable.h"

std::vector<std::thread> procVec;
std::mutex threadMutex;

void rem_unfix_pe(SRA sra, long long int ram_b) {
  std::string inFile1Str(sra.get_sra_path_corr().first.c_str());
  std::string inFile2Str(sra.get_sra_path_corr().second.c_str());
  std::string outFile1Str(std::string(sra.get_sra_path_corr().first.replace_extension("fixB.fq").c_str()));
  std::string outFile2Str(std::string(sra.get_sra_path_corr().second.replace_extension("fixB.fq").c_str()));
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

  long long int ram_b_per_file = ram_b / 2;

  inFile1.seekg(0, inFile1.end);
  long int lenFile1 = inFile1.tellg();
  inFile1.seekg(0, inFile1.beg);

  inFile2.seekg(0, inFile2.end);
  long int lenFile2 = inFile2.tellg();
  inFile2.seekg(0, inFile2.beg);

  char * inFile1Data = new char[ram_b_per_file];
  char * inFile2Data = new char[ram_b_per_file];
  std::streamsize s1;
  std::streamsize s2;

  char * currPos1;
  char * currPos2;
  char * nlPos1;
  char * nlPos2;
  char * nlPos1Prev;
  char * nlPos2Prev;
  char * writeStart1;
  char * writeStart2;
  char * writeEnd1;
  char * writeEnd2;
  char * inFile1L;
  char * inFile2L;

  std::string readMatch;
  std::string read;
  while (!inFile1.eof() && !inFile2.eof()) {
    inFile1.read(inFile1Data, ram_b_per_file);
    inFile2.read(inFile2Data, ram_b_per_file);

    s1 = inFile1.gcount();
    s2 = inFile2.gcount();

    currPos1 = inFile1Data;
    currPos2 = inFile2Data;
    nlPos1 = inFile1Data;
    nlPos2 = inFile2Data;
    writeStart1 = inFile1Data;
    writeStart2 = inFile2Data;

    inFile1L = inFile1Data + s1;
    inFile2L = inFile2Data + s2;

    // Unget character until end of inFile1 / inFile2 buffer before a read
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

    // inFile1 buffer in correct position
    // Correctly position the end of inFile2 buffer
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
      if (read > readMatch) {
        while (read.compare(readMatch)) {
          // Unget file2 until '@' or '>' reached
          // Decrement s2 accordingly
          // Remove last char from inFile2Data accordingly
          // Get new read string
        }
      }
      else if (read < readMatch) {
        while (read.compare(readMatch)) {
          // Get file2 until '@' or '>' reached
          // Increment s2 accordingly
          // Add last char to inFile2Data accordingly
          // Get new read string
        }
      }
      else {
        // Buffer in position -- do nothing
      }
    }

    // inFile2 buffer in correct position
    // Correctly position the end of inFile1 buffer
    else {
      inFile2.get();
      inFile2 >> readMatch;
      while (inFile2.peek() != '@' && inFile2.peek() != '>') {
        inFile2.unget();
      }
      while (inFile1.peek() != '@' && inFile2.peek() != '>') {
        inFile1.unget();
        inFile1Data[s1 - 1] = '\0';
        s1--;
      }
      inFile1.get();
      inFile1 >> read;
      if (read > readMatch) {
        while (read.compare(readMatch)) {
          // Unget file1 until '@' or '>' reached
          // Decrement s1 accordingly
          // Remove last char from inFile1Data accordingly
          // Get new read string
        }
      }
      else if (read < readMatch) {
        while (read.compare(readMatch)) {
          // Get file1 until '@' or '>' reached
          // Increment s1 accordingly
          // Add last char to inFile2Data accordingly
          // Get new read string
        }
      }
      else {
        // Buffer in position -- do nothing
      }
    }
    
    while (nlPos1 != inFile1L && nlPos2 != inFile2L) {
      nlPos1Prev = nlPos1;
      nlPos2Prev = nlPos2;
      nlPos1 = std::find(nlPos1 + 1, inFile1L, '\n');
      nlPos2 = std::find(nlPos2 + 1, inFile2L, '\n');
      if (strncmp(nlPos1 - 5, "error", 5) == 0 ||
          strncmp(nlPos2 - 5, "error", 5) == 0) {
        writeEnd1 = nlPos1Prev;
        outFile1.write(writeStart1, writeEnd1 - writeStart1);
        writeEnd2 = nlPos2Prev;
        outFile2.write(writeStart2, writeEnd2 - writeStart2);

        for (int i = 0; i < 3; i++) {
          nlPos1 = std::find(nlPos1 + 1, inFile1L, '\n');
          nlPos2 = std::find(nlPos2 + 1, inFile2L, '\n');
        }
        writeStart1 = nlPos1;
        writeStart2 = nlPos2;
      }
    }
    outFile1.write(writeStart1, inFile1Data + s1 - writeStart1);
    outFile2.write(writeStart2, inFile2Data + s2 - writeStart2);
  }
/*
  std::string currLine1;
  std::string currLine1;
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
*/
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
