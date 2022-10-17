#include "rem_overrep.h"

std::pair<std::vector<std::string>, std::vector<std::string>> get_overrep_seqs_pe(SRA sra) {
  std::vector<std::string> overrepSeqs1;
  std::vector<std::string> overrepSeqs2;

  std::string inFile1Str(std::string(sra.get_fastqc_dir().first.c_str()) + "_fastqc.html");
  std::string inFile2Str(std::string(sra.get_fastqc_dir().second.c_str()) + "_fastqc.html");

  std::ifstream inFile1(inFile1Str);
  std::ifstream inFile2(inFile2Str);

  boost::regex rgx("(?<=<tr><td>)[ATCG]*(?=\<\/td>)");
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
    overrepSeqs1.push_back(res.str());
    sStart = res.suffix().first;
  }

  sStart = inFile2Data.begin();
  sEnd = inFile2Data.end();

  while (boost::regex_search(sStart, sEnd, res, rgx)) {
    overrepSeqs2.push_back(res.str());
    sStart = res.suffix().first;
  }

  std::pair<std::vector<std::string>, std::vector<std::string>> overrepSeqs(overrepSeqs1, overrepSeqs2);  

  return overrepSeqs;
}


std::vector<std::string> get_overrep_seqs_se(SRA sra) {
  std::vector<std::string> overrepSeqs;

  std::string inFileStr(std::string(sra.get_fastqc_dir().first.c_str()) + "_fastqc.html");

  std::ifstream inFile(inFileStr);

  boost::regex rgx("(?<=<tr><td>)[ATCG]*(?=\<\/td>)");
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

