#include "rem_unfixable.h"


void rem_unfix_pe(SRA sra, long long int ram_b) {
  std::string inFile1Str(sra.get_sra_path_corr().first.c_str());
  std::string inFile2Str(sra.get_sra_path_corr().second.c_str());
  std::string outFile1Str(std::string(sra.get_sra_path_corr().first.replace_extension("fix.fq").c_str()));
  std::string outFile2Str(std::string(sra.get_sra_path_corr().second.replace_extension("fix.fq").c_str()));
  std::ifstream inFile1(inFile1Str);
  std::ifstream inFile2(inFile2Str);
  std::ofstream outFile1(outFile1Str);
  std::ofstream outFile2(outFile2Str);

  long long int ram_b_per_file = ram_b / 2;

  std::string inFile1Data;
  std::string inFile2Data;

  inFile1Data.reserve(ram_b_per_file);
  inFile2Data.reserve(ram_b_per_file);

  std::streamsize s1;
  std::streamsize s2;

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

  while (!inFile1.eof() || !inFile2.eof()) {
    inFile1.read(&inFile1Data[0], ram_b_per_file);
    inFile2.read(&inFile2Data[0], ram_b_per_file);

    s1 = inFile1.gcount();
    s2 = inFile2.gcount();

    nlPos1 = &inFile1Data[0];
    nlPos2 = &inFile2Data[0];
    writeStart1 = &inFile1Data[0];
    writeStart2 = &inFile2Data[0];

    inFile1L = &inFile1Data[0] + s1;
    inFile2L = &inFile2Data[0] + s2;

    align_file_buffer(inFile1, inFile2, &inFile1Data[0], &inFile2Data[0], s1, s2);

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
    outFile1.write(writeStart1, &inFile1Data[0] + s1 - writeStart1);
    outFile2.write(writeStart2, &inFile2Data[0] + s2 - writeStart2);
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

  inFile.seekg(0, inFile.end);
  long int lenFile = inFile.tellg();
  inFile.seekg(0, inFile.beg);

  char * inFileData = new char[ram_b];
  std::streamsize s;

  char * nlPos;
  char * nlPosPrev;
  char * writeStart;
  char * writeEnd;
  char * inFileL;

  while (!inFile.eof()) {
    inFile.read(inFileData, ram_b);
    s = inFile.gcount();
    nlPos = inFileData;
    writeStart = inFileData;
    inFileL = inFileData + s;

    if (!inFile.eof()) {
      while (inFile.peek() != '@' && inFile.peek() != '>') {
        inFile.unget();
        inFileData[s - 1] = '\0';
        s--;
      }
    }

    while (nlPos != inFileL) {
      nlPosPrev = nlPos;
      nlPos = std::find(nlPos + 1, inFileL, '\n');
      if (strncmp(nlPos - 5, "error", 4) == 0) {
        writeEnd = nlPosPrev;
        outFile.write(writeStart, writeEnd - writeStart);

        for (int i = 0; i < 3; i++) {
          nlPos = std::find(nlPos + 1, inFileL, '\n');
        }
        writeStart = nlPos;
      }
    }
    outFile.write(writeStart, inFileData + s - writeStart);
  }
  inFile.close();
  outFile.close();
}


void rem_unfix_bulk(const std::vector<SRA> & sras, std::string ram_gb, std::string logFile) {
  std::cout << "\nRemoving unfixable reads for:\n" << std::endl;
  summarize_all_sras(sras, logFile, 2);
  long long int ram_b = (long long int)stoi(ram_gb) * 1000000000;
  for (auto sra : sras) {
    if (fs::exists(sra.get_sra_path_corr_fix().first.c_str())) {
      std::cout << "Fixed version found for: " << sra.get_accession() << std::endl;
      continue;
    }
    if (sra.is_paired()) {
      rem_unfix_pe(sra, ram_b);
    }
    else {
      rem_unfix_se(sra, ram_b);
    }
  }
}
