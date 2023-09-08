#include "rem_overrep.h"

// Given a paired-end SRA run object post-FastQC, return a pair of vectors containing its
// overrepresented sequences
std::pair<std::vector<std::string>, std::vector<std::string>> get_overrep_seqs_pe(SRA sra) {

  std::pair<std::vector<std::string>, std::vector<std::string>> overrepSeqs;

  std::string inFile1Str(std::string(sra.get_fastqc_dir_2().first.c_str()) + ".filt_fastqc.html");
  std::string inFile2Str(std::string(sra.get_fastqc_dir_2().second.c_str()) + ".filt_fastqc.html");

  std::ifstream inFile1(inFile1Str);
  std::ifstream inFile2(inFile2Str);

  boost::regex rgx("(?<=<tr><td>)[ATCG]*(?=\\<\\/td>)");
  boost::smatch res;

  std::string inFile1Data;
  while (getline(inFile1, inFile1Data));
  std::string inFile2Data;
  while (getline(inFile2, inFile2Data));

  inFile1.close();
  inFile2.close();

  std::string::const_iterator sStart, sEnd;

  sStart = inFile1Data.begin();
  sEnd = inFile1Data.end();

  while (boost::regex_search(sStart, sEnd, res, rgx)) {
    overrepSeqs.first.push_back(res.str());
    sStart = res.suffix().first;
  }

  sStart = inFile2Data.begin();
  sEnd = inFile2Data.end();

  while (boost::regex_search(sStart, sEnd, res, rgx)) {
    overrepSeqs.second.push_back(res.str());
    sStart = res.suffix().first;
  }

  return overrepSeqs;
}

// Given a single-end SRA run object post-FastQC, return a vector containing its overrepresented
// sequences
std::vector<std::string> get_overrep_seqs_se(SRA sra) {
  std::vector<std::string> overrepSeqs;

  std::string inFileStr(std::string(sra.get_fastqc_dir_2().first.c_str()) + ".filt_fastqc.html");

  std::ifstream inFile(inFileStr);

  boost::regex rgx("(?<=<tr><td>)[ATCG]*(?=\\<\\/td>)");
  boost::smatch res;

  std::string inFileData;
  while (getline(inFile, inFileData));

  inFile.close();

  std::string::const_iterator sStart, sEnd;

  sStart = inFileData.begin();
  sEnd = inFileData.end();

  while (boost::regex_search(sStart, sEnd, res, rgx)) {
    overrepSeqs.push_back(res.str());
    sStart = res.suffix().first;
  }

  return overrepSeqs;
}

// Given a paired-end SRA run's sequence data files and vectors containing its
// overrepresented sequences, parse the files and remove reads containing these sequences
bool rem_overrep_pe(std::pair<std::string, std::string> sraRunIn,
                    std::pair<std::string, std::string> sraRunOut,
                    uintmax_t ram_b, bool dispOutput, bool compressFiles,
                    std::pair<std::vector<std::string>, std::vector<std::string>> overrepSeqs,
                    std::string logFile) {
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

  if (overrepSeqs.first.empty() && overrepSeqs.second.empty()) {
    system(std::string("cp " + sraRunIn.first + " " + sraRunOut.first).c_str());
    system(std::string("cp " + sraRunIn.second + " " + sraRunOut.second).c_str());
    return true;
  }

  size_t lenOverrep1 = 0;
  size_t lenOverrep2 = 0;

  if (!overrepSeqs.first.empty()) {
    lenOverrep1 = overrepSeqs.first.front().length();
  }
  if (!overrepSeqs.second.empty()) {
    lenOverrep2 = overrepSeqs.second.front().length();
  }

  uintmax_t ram_b_per_file = ram_b / 2;

  //std::string inFile1Data;
  //std::string inFile2Data;

  //inFile1Data.reserve(ram_b_per_file);
  //inFile2Data.reserve(ram_b_per_file);
  char * inFile1Data;
  char * inFile2Data;

  inFile1Data = new char[ram_b_per_file];
  inFile2Data = new char[ram_b_per_file];

  std::streamsize s1;
  std::streamsize s2;

  char * nlPos1;
  char * nlPos2;
  char * nlPos1Head;
  char * nlPos2Head;
  char * nlPos1Prev;
  char * nlPos2Prev;
  char * writeStart1;
  char * writeStart2;
  char * writeEnd1;
  char * writeEnd2;
  char * inFile1L;
  char * inFile2L;
  int seqLength;
  bool overRep;
 
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

  uintmax_t numOrep = 0;
  uintmax_t numReads = 0;
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
      nlPos1Head = std::find(nlPos1 + 1, inFile1L, '\n');
      nlPos2Head = std::find(nlPos2 + 1, inFile2L, '\n');
      nlPos1 = std::find(nlPos1Head + 1, inFile1L, '\n');
      nlPos2 = std::find(nlPos2Head + 1, inFile2L, '\n');

      seqLength = nlPos1 - nlPos1Head - 1;
      overRep = false;
      if (!overrepSeqs.first.empty()) {
        for (int i = 0; i < seqLength - lenOverrep1; i++) {
          for (auto const & seq : overrepSeqs.first) {
            if (strncmp(nlPos1Head + 1 + i, seq.c_str(), lenOverrep1) == 0) {
              overRep = true;
              goto checkOverrep;
            }
          }
        }
      }
      if (!overrepSeqs.second.empty()) {
        for (int i = 0; i < seqLength - lenOverrep2; i++) {
          for (auto const & seq : overrepSeqs.second) {
            if (strncmp(nlPos2Head + 1 + i, seq.c_str(), lenOverrep2) == 0) {
              overRep = true;
              goto checkOverrep;
            }
          }
        }
      }
      checkOverrep:if (overRep) {
        numOrep++;
        readName = std::string(nlPos1Prev + 1, nlPos1Head);
        if (dispOutput) {
          logOutput("Overrepresented read removed: " + readName.substr(1, readName.find(" ")), logFile);
        }
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
 
        numReads++;
        for (int i = 0; i < 2; i++) {
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
        numReads++;
        okReads++;
        for (int i = 0; i < 2; i++) {
          nlPos1 = std::find(nlPos1 + 1, inFile1L, '\n');
          nlPos2 = std::find(nlPos2 + 1, inFile2L, '\n');
        }
      }
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
    return false;
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

  if (dispOutput) {
    std::stringstream percentStream;
    percentStream << std::fixed << std::setprecision(2)
                  << (float(numOrep) * 100) / float(numOrep + numReads);
    logOutput("Removal of overrepresented reads finished", logFile);
    logOutput("\n  SUMMARY", logFile);
    logOutput("    " + std::to_string(numOrep) + " reads removed (" +
              percentStream.str() + " %)", logFile);
    percentStream.str("");
    percentStream << std::fixed << std::setprecision(2)
                  << (float(numReads) * 100) / (float(numOrep + numReads));
    logOutput("    " + std::to_string(numReads - numOrep) + " reads retained (" +
              percentStream.str() + " %)", logFile);
  }
  return true;
}

// Given a single-end SRA run's sequence data file and a vector containing its
// overrepresented sequences, parse the file and remove reads containing these sequences
bool rem_overrep_se(std::string sraRunIn, std::string sraRunOut,
                    uintmax_t ram_b, bool dispOutput, bool compressFiles,
                    std::vector<std::string> overrepSeqs,
                    std::string logFile) {
  std::ifstream inFile(sraRunIn);
  std::ofstream outFile(sraRunOut);

  if (overrepSeqs.empty()) {
    system(std::string("cp " + sraRunIn + " " + sraRunOut).c_str());
    return true;
  }

  size_t lenOverrep = overrepSeqs.front().size();

  char * inFileData;

  inFileData = new char[ram_b];

  std::streamsize s;

  char * nlPos;
  char * nlPosHead;
  char * nlPosPrev;
  char * writeStart;
  char * writeEnd;
  char * inFileL;
  int seqLength;
  bool overRep;

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

  uintmax_t numOrep = 0;
  uintmax_t numReads = 0;
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

    uintmax_t okReads = 0;
    while (nlPos != inFileL) {
      nlPosPrev = nlPos;
      nlPosHead = std::find(nlPos + 1, inFileL, '\n');
      nlPos = std::find(nlPosHead + 1, inFileL, '\n');
      seqLength = nlPos - nlPosHead - 1;
      overRep = false;
      for (int i = 0; i < seqLength - lenOverrep; i++) {
        for (auto const & seq : overrepSeqs) {
          if (strncmp(nlPosHead + 1 + i, &seq[0], lenOverrep) == 0) {
            overRep = true;
            goto checkOverrep;
          }
        }
      }
      checkOverrep:if (overRep) {
        numOrep++;
        readName = std::string(nlPosPrev + 1, nlPosHead);
        if (dispOutput) {
          logOutput("Overrepresented read removed: " + readName.substr(1, readName.find(" ")), logFile);
        }
        if (okReads == 0) {
          writeEnd = nlPosPrev + 1;
        }
        else {
          writeEnd = nlPosPrev;
        }
        if (compressFiles) {
          outputStream.write(writeStart, writeEnd - writeStart);
        }
        else {
          outFile.write(writeStart, writeEnd - writeStart);
        }

        numReads++;
        for (int i = 0; i < 2; i++) {
          nlPos = std::find(nlPos + 1, inFileL, '\n');
        }
        if (okReads == 0) {
          writeStart = nlPos + 1;
        }
        else {
          writeStart = nlPos;
        }
      }
      else {
        numReads++;
        okReads++;
        for (int i = 0; i < 2; i++) {
          nlPos = std::find(nlPos + 1, inFileL, '\n');
        }
      }
    }
    if (compressFiles) {
      outputStream.write(writeStart, writeEnd - writeStart);
    }
    else {
      outFile.write(writeStart, inFileData + s - writeStart);
    }
  }
  if (!outFile.good() || !outputStream.good()) {
    return false;
  }
  if (!compressFiles) {
    inFile.close();
  }
  else {
    io::close(gzInBuffer);
    io::close(gzOutBuffer);
  }

  delete [] inFileData;
  outFile.close();

  if (dispOutput) {
    std::stringstream percentStream;
    percentStream << std::fixed << std::setprecision(2)
                  << (float(numOrep) * 100) / float(numOrep + numReads);
    logOutput("Removal of overrepresented reads finished", logFile);
    logOutput("\n  SUMMARY", logFile);
    logOutput("    " + std::to_string(numOrep) + " reads removed (" +
              percentStream.str() + " %)", logFile);
    percentStream.str("");
    percentStream << std::fixed << std::setprecision(2)
                  << (float(numReads) * 100) / (float(numOrep + numReads));
    logOutput("    " + std::to_string(numReads - numOrep) + " reads retained (" +
              percentStream.str() + " %)", logFile);
  }
  return true;
}

