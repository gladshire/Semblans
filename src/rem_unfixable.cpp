
#include "rem_unfixable.h"

// Given a paired-end SRA run's sequence data files post-Rcorrector,
// remove all reads Rcorrector flagged as "unfixable error"
long int rem_unfix_pe(std::pair<std::string, std::string> sraRunIn,
                      std::pair<std::string, std::string> sraRunOut,
                      uintmax_t ram_b, bool dispOutput, bool compressFiles,
                      std::string logFile) {
  std::ifstream inFile1;
  std::ifstream inFile2;
  std::ofstream outFile1(sraRunOut.first);
  std::ofstream outFile2(sraRunOut.second);
  std::stringstream percentStream;

  if (compressFiles) {
    inFile1.open(sraRunIn.first, std::ios_base::binary);
    inFile2.open(sraRunIn.second, std::ios_base::binary);
  }
  else {
    inFile1.open(sraRunIn.first);
    inFile2.open(sraRunIn.second);
  }

  long long int ram_b_per_file = ram_b / 2;

  char * inFile1Data;
  char * inFile2Data;

  inFile1Data = new char[ram_b_per_file];
  inFile2Data = new char[ram_b_per_file];

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

  uintmax_t numUnfix = 0;
  uintmax_t numUnfix100k = 0;
  uintmax_t numReads = 0;
  float percentUnfix;
  std::string readName;
  while  ((!inputStream1.eof() && !inputStream2.eof() && !inFile1.eof() && !inFile2.eof()) &&
          (inputStream1.good() && inputStream2.good() && inFile1.good() && inFile2.good())) {

    if (compressFiles) {
      inputStream1.read(inFile1Data, ram_b_per_file);
      inputStream2.read(inFile2Data, ram_b_per_file);

      s1 = inputStream1.gcount();
      s2 = inputStream2.gcount();
    }
    else {
      inFile1.read(inFile1Data, ram_b_per_file);
      inFile2.read(inFile2Data, ram_b_per_file);

      s1 = inFile1.gcount();
      s2 = inFile2.gcount();
    }

    nlPos1 = inFile1Data;
    nlPos2 = inFile2Data;
    writeStart1 = inFile1Data;
    writeStart2 = inFile2Data;

    inFile1L = inFile1Data + s1;
    inFile2L = inFile2Data + s2;

    if (compressFiles) {
      align_file_buffer(inputStream1, inputStream2, inFile1Data, inFile2Data, s1, s2);
    }
    else {
      align_file_buffer(inFile1, inFile2, inFile1Data, inFile2Data, s1, s2);
    }

    uintmax_t okReads = 0;
    while (nlPos1 != inFile1L && nlPos2 != inFile2L) {
      nlPos1Prev = nlPos1;
      nlPos2Prev = nlPos2;
      nlPos1 = std::find(nlPos1 + 1, inFile1L, '\n');
      nlPos2 = std::find(nlPos2 + 1, inFile2L, '\n');
      if (strncmp(nlPos1 - 5, "error", 5) == 0 ||
          strncmp(nlPos2 - 5, "error", 5) == 0) {
        numUnfix++;
        if (okReads == 0) {
          writeEnd1 = nlPos1Prev + 1;
          writeEnd2 = nlPos2Prev + 1;
        }
        else {
          writeEnd1 = nlPos1Prev;
          writeEnd2 = nlPos2Prev;
        }

        if (okReads > 0) {
          if (compressFiles) {
            outputStream1.write(writeStart1, writeEnd1 - writeStart1);
            outputStream2.write(writeStart2, writeEnd2 - writeStart2);
          }
          else {
            outFile1.write(writeStart1, writeEnd1 - writeStart1);
            outFile2.write(writeStart2, writeEnd2 - writeStart2);
          }
        }
        if (numUnfix % 100000 == 0) {
          numUnfix100k++;
          if (dispOutput) {
            std::cout << "\r  Unfixable reads removed: " + std::to_string(numUnfix100k * 100000) +
                      " ...     " << std::flush;
          }
        }

        for (int i = 0; i < 3; i++) {
          nlPos1 = std::find(nlPos1 + 1, inFile1L, '\n');
          nlPos2 = std::find(nlPos2 + 1, inFile2L, '\n');
        }
        if (okReads == 0) {
          writeStart1 = nlPos1 + 1;
          writeStart2 = nlPos2 + 1;
        }
        else {
          writeStart1 = nlPos1;
          writeStart2 = nlPos2;
        }
      }
      else {
        okReads++;
        for (int i = 0; i < 3; i++) {
          nlPos1 = std::find(nlPos1 + 1, inFile1L, '\n');
          nlPos2 = std::find(nlPos2 + 1, inFile2L, '\n');
        }
      }
      numReads++;
    }
    if (compressFiles) {
      outputStream1.write(writeStart1, inFile1Data + s1 - writeStart1);
      outputStream2.write(writeStart2, inFile2Data + s2 - writeStart2);
    }
    else {
      outFile1.write(writeStart1, inFile1Data + s1 - writeStart1);
      outFile2.write(writeStart2, inFile2Data + s2 - writeStart2);
    }
  }
  if (!outputStream1.good() || !outputStream2.good() || !outFile1.good() || !outFile2.good()) {
    std::cerr << "ERROR: Writing output failed for:\n  "
              << sraRunIn.first << "\n  " << sraRunIn.second << std::endl;
    return -1;
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
  delete [] inFile1Data;
  delete [] inFile2Data;
  outFile1.close();
  outFile2.close();

  std::string percentStr;
  if (dispOutput) {
    percentUnfix = (float(numUnfix) * 100) / float(numReads);
    logOutput("\nRemoval of unfixable reads finished\n", logFile);
    logOutput("\n  SUMMARY", logFile);
    logOutput("\n    Reads removed: " + std::to_string(numUnfix) +
              "(" + getPercent(percentUnfix, 2) + "%)", logFile);
    logOutput("\n    Reads retained: " + std::to_string(numReads - numUnfix) +
              "(" + getPercent(100.0 - percentUnfix, 2) + "%)\n", logFile);
  }
  return numUnfix;
}

// Given a single-end SRA run's sequence data file post-Rcorrector,
// remove all reads Rcorrector flagged as "unfixable error"
long int rem_unfix_se(std::string sraRunIn, std::string sraRunOut,
                      uintmax_t ram_b, bool dispOutput, bool compressFiles,
                      std::string logFile) {

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

  uintmax_t numUnfix = 0;
  uintmax_t numUnfix100k = 0;
  uintmax_t numReads = 0;
  float percentUnfix;
  std::string readName;
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

    uintmax_t okReads = 0;
    while (nlPos != inFileL) {
      nlPosPrev = nlPos;
      nlPos = std::find(nlPos + 1, inFileL, '\n');
      if (strncmp(nlPos - 5, "error", 4) == 0) {
        numUnfix++;
        writeEnd = nlPosPrev;
        if (okReads > 0) {
          if (compressFiles) {
            outputStream.write(writeStart, writeEnd - writeStart);
          }
          else {
            outFile.write(writeStart, writeEnd - writeStart);
          }
        }
        if (numUnfix % 100000 == 0) {
          numUnfix100k++;
          std::cout << "\r  Unfixable reads removed: " + std::to_string(numUnfix100k * 100000) +
                    " ...     " << std::flush;
        }
        for (int i = 0; i < 3; i++) {
          nlPos = std::find(nlPos + 1, inFileL, '\n');
        }
        writeStart = nlPos;
      }
      else {
        okReads++;
        for (int i = 0; i < 3; i++) {
          nlPos = std::find(nlPos + 1, inFileL, '\n');
        }
      }
      numReads++;
    }
    if (!outFile.good() || !outputStream.good()) {
      std::cerr << "ERROR: Writing output failed for:\n  " << sraRunIn << std::endl;
      return -1;
    }
    logOutput("\r  Unfixable reads removed: " + std::to_string(numUnfix) +
              "         ", logFile);
    if (compressFiles) {
      outputStream.write(writeStart, writeEnd + s - writeStart);
    }
    else {
      outFile.write(writeStart, inFileData + s - writeStart);
    }
  }

  delete [] inFileData;

  if (!compressFiles) {
    inFile.close();
  }
  else {
    io::close(gzInBuffer);
    io::close(gzOutBuffer);
  }
  outFile.close();

  if (dispOutput) {
    percentUnfix = (float(numUnfix) * 100) / float(numReads);
    logOutput("\nRemoval of unfixable reads finished", logFile);
    logOutput("\n  SUMMARY", logFile);
    logOutput("\n    " + std::to_string(numUnfix) + " reads removed (" +
              getPercent(percentUnfix, 2) + "%)", logFile);
    logOutput("\n    " + std::to_string(numReads - numUnfix) + " reads retained (" +
              getPercent(100.0 - percentUnfix, 2) + "%)\n", logFile);
  }
  return numUnfix;
}
