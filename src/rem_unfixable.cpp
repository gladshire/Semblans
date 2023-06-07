#include "rem_unfixable.h"


void rem_unfix_pe(std::pair<std::string, std::string> sraRunIn,
                  std::pair<std::string, std::string> sraRunOut,
                  uintmax_t ram_b, bool compressFiles) {
  std::ifstream inFile1;
  std::ifstream inFile2;
  std::ofstream outFile1(sraRunOut.first);
  std::ofstream outFile2(sraRunOut.second);

  if (compressFiles) {
    inFile1.open(sraRunIn.first, std::ios_base::binary);
    inFile2.open(sraRunIn.second, std::ios_base::binary);
  }
  else {
    inFile1.open(sraRunIn.first);
    inFile2.open(sraRunIn.second);
  }

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

  io::filtering_streambuf<io::input> gzInBuffer1;
  io::filtering_streambuf<io::input> gzInBuffer2;
   
  gzInBuffer1.push(io::gzip_decompressor());
  gzInBuffer2.push(io::gzip_decompressor());
    
  gzInBuffer1.push(inFile1);
  gzInBuffer2.push(inFile2);

  std::istream inputStream1(&gzInBuffer1);
  std::istream inputStream2(&gzInBuffer2);

  while (!inFile1.eof() || !inFile2.eof()) {
    if (compressFiles) {
      inputStream1.read(&inFile1Data[0], ram_b_per_file);
      inputStream2.read(&inFile2Data[0], ram_b_per_file);
  
      s1 = inputStream1.gcount();
      s2 = inputStream2.gcount();
    }
    else {
      inFile1.read(&inFile1Data[0], ram_b_per_file);
      inFile2.read(&inFile2Data[0], ram_b_per_file);
    
      s1 = inFile1.gcount();
      s2 = inFile2.gcount();
    }

    nlPos1 = &inFile1Data[0];
    nlPos2 = &inFile2Data[0];
    writeStart1 = &inFile1Data[0];
    writeStart2 = &inFile2Data[0];

    inFile1L = &inFile1Data[0] + s1;
    inFile2L = &inFile2Data[0] + s2;
    
    if (compressFiles) {
      align_file_buffer(inputStream1, inputStream2, &inFile1Data[0], &inFile2Data[0], s1, s2);
    }
    else {
      align_file_buffer(inFile1, inFile2, &inFile1Data[0], &inFile2Data[0], s1, s2);
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
    outFile1.write(writeStart1, &inFile1Data[0] + s1 - writeStart1);
    outFile2.write(writeStart2, &inFile2Data[0] + s2 - writeStart2);
  }
  if (!compressFiles) {
    inFile1.close();
    inFile2.close();
  }
  outFile1.close();
  outFile2.close();
}


void rem_unfix_se(std::string sraRunIn, std::string sraRunOut,
                  uintmax_t ram_b) {

  std::ifstream inFile(sraRunIn);
  std::ofstream outFile(sraRunOut);

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
