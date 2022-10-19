#include "rem_unfixable.h"

std::vector<std::thread> procVec;
std::mutex threadMutex;


// Notes:
//  - Concurrent file writing bad
//  - Fastest file reading/writing:
//    - Use buffer in RAM
//    - Find optimal size, based on user's choice
//      - char * limit: 65535
//      - Possibly use strings?
//        - Would need to adapt algo to work on strings

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

  inFile1.seekg(0, inFile1.end);
  long int lenFile1 = inFile1.tellg();
  inFile1.seekg(0);

  inFile2.seekg(0, inFile2.end);
  long int lenFile2 = inFile2.tellg();
  inFile2.seekg(0);

  char * inFile1Data;
  char * inFile2Data;

  long long int writeSize;

  inFile1Data = new char[ram_b_per_file];
  inFile2Data = new char[ram_b_per_file];
  writeSize = ram_b_per_file;


  /*int writeSize = 4096;

  char inFile1Data[writeSize];
  char inFile2Data[writeSize];
*/
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
  while (!inFile1.eof() || !inFile2.eof()) {
    inFile1.read(inFile1Data, writeSize);
    inFile2.read(inFile2Data, writeSize);

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
          inFile2.get();
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

  char * currPos;
  char * nlPos;
  char * nlPosPrev;
  char * writeStart;
  char * writeEnd;
  char * inFileL;

  while (!inFile.eof()) {
    inFile.read(inFileData, ram_b);
    s = inFile.gcount();
    currPos = inFileData;
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


void rem_unfix_bulk(std::vector<SRA> sras, std::string ram_gb) {
  std::cout << "\nRemoving unfixable reads for:\n" << std::endl;
  summarize_all_sras(sras);
  long long int ram_b = stoi(ram_gb) * 1000000000;
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
