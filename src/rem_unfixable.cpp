// TODO: Not compressing files properly. Missing last few lines

#include "rem_unfixable.h"

// Given a paired-end SRA run's sequence data files post-Rcorrector,
// remove all reads Rcorrector flagged as "unfixable error"
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

  // Instantiate input stream for compressed file
  io::filtering_streambuf<io::input> gzInBuffer1;
  io::filtering_streambuf<io::input> gzInBuffer2;
     
  gzInBuffer1.push(io::gzip_decompressor());
  gzInBuffer2.push(io::gzip_decompressor());
    
  gzInBuffer1.push(inFile1);
  gzInBuffer2.push(inFile2);

  std::istream inputStream1(&gzInBuffer1);
  std::istream inputStream2(&gzInBuffer2);


  // Instantiate output stream for compressed file
  io::filtering_streambuf<io::output> gzOutBuffer1;
  io::filtering_streambuf<io::output> gzOutBuffer2;

  gzOutBuffer1.push(io::gzip_compressor(io::gzip_params(io::gzip::best_speed)));
  gzOutBuffer2.push(io::gzip_compressor(io::gzip_params(io::gzip::best_speed)));

  gzOutBuffer1.push(outFile1);
  gzOutBuffer2.push(outFile2);

  std::ostream outputStream1(&gzOutBuffer1);
  std::ostream outputStream2(&gzOutBuffer2);
 
  int numUnfix = 0;
 
  while  ((!inputStream1.eof() && !inputStream2.eof() && !inFile1.eof() && !inFile2.eof()) &&
         (inputStream1.good() && inputStream2.good() && inFile1.good() && inFile2.good())) {

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
        numUnfix++;
        if (numUnfix == 1) {
          writeEnd1 = nlPos1Prev - 1;
          writeEnd2 = nlPos2Prev - 1;
        }
        else {
          writeEnd1 = nlPos1Prev;
          writeEnd2 = nlPos2Prev;
        }
 
        if (compressFiles) {
          outputStream1.write(writeStart1, writeEnd1 - writeStart1);
          outputStream2.write(writeStart2, writeEnd2 - writeStart2);
        }
        else {
          outFile1.write(writeStart1, writeEnd1 - writeStart1);
          outFile2.write(writeStart2, writeEnd2 - writeStart2);
        }

        for (int i = 0; i < 3; i++) {
          nlPos1 = std::find(nlPos1 + 1, inFile1L, '\n');
          nlPos2 = std::find(nlPos2 + 1, inFile2L, '\n');
        }
        writeStart1 = nlPos1;
        writeStart2 = nlPos2;
      }
    }
    if (compressFiles) {
      outputStream1.write(writeStart1, &inFile1Data[0] + s1 - writeStart1);
      outputStream2.write(writeStart2, &inFile2Data[0] + s2 - writeStart2);
    }
    else {
      outFile1.write(writeStart1, &inFile1Data[0] + s1 - writeStart1);
      outFile2.write(writeStart2, &inFile2Data[0] + s2 - writeStart2);
    }
  }
  if (!compressFiles) {
    inFile1.close();
    inFile2.close();
  }
  else {
    io::close(gzInBuffer1);
    io::close(gzInBuffer2);
    io::close(gzOutBuffer1);
    io::close(gzOutBuffer2);
  }
  outFile1.close();
  outFile2.close();
}

// Given a single-end SRA run's sequence data file post-Rcorrector,
// remove all reads Rcorrector flagged as "unfixable error"
void rem_unfix_se(std::string sraRunIn, std::string sraRunOut,
                  uintmax_t ram_b, bool compressFiles) {

  std::ifstream inFile;
  std::ofstream outFile(sraRunOut);

  if (compressFiles) {
    inFile.open(sraRunIn, std::ios_base::binary);
  }
  else {
    inFile.open(sraRunIn);
  }

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

  // Instantiate input stream of compressed file
  io::filtering_streambuf<io::input> gzInBuffer;
  gzInBuffer.push(io::gzip_decompressor());
  gzInBuffer.push(inFile);
  std::istream inputStream(&gzInBuffer);

  // Instantiate output stream to compressed file
  io::filtering_streambuf<io::output> gzOutBuffer;
  gzOutBuffer.push(io::gzip_compressor(io::gzip_params(io::gzip::best_speed)));
  gzOutBuffer.push(outFile);
  std::ostream outputStream(&gzOutBuffer);

  while ((!inputStream.eof() && !inFile.eof()) && (inputStream.good() && inFile.good())) {
    if (compressFiles) {
      inputStream.read(inFileData, ram_b);
      s = inputStream.gcount();
    }
    else {
      inFile.read(inFileData, ram_b);
      s = inFile.gcount();
    }

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
        if (compressFiles) {
          outputStream.write(writeStart, writeEnd - writeStart);
        }
        else {
          outFile.write(writeStart, writeEnd - writeStart);
        }

        for (int i = 0; i < 3; i++) {
          nlPos = std::find(nlPos + 1, inFileL, '\n');
        }
        writeStart = nlPos;
      }
    }
    if (compressFiles) {
      outputStream.write(writeStart, writeEnd + s - writeStart);
    }
    else {
      outFile.write(writeStart, inFileData + s - writeStart);
    }
  }
  if (!compressFiles) {
    inFile.close();
  }
  else {
    io::close(gzInBuffer);
    io::close(gzOutBuffer);
  }
  outFile.close();
}
