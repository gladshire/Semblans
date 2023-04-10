#include "rem_chimera.h"

std::vector<std::string> getLineVec(std::string line) {
  std::vector<std::string> lineVec;
  
  if (line.length() < 3) {
    return lineVec;
  }
  std::string delim = "\t";
  std::string ttok;
  size_t tpos = 0;
  while ((tpos = line.find(delim)) != std::string::npos) {
    ttok = line.substr(0, tpos);
    lineVec.push_back(ttok);
    line.erase(0, tpos + 1);
  }
  lineVec.push_back(line);
  return lineVec;
}


double getQcov(std::vector<std::string> currLineVec) {
  double start = std::stod(currLineVec[10]);
  double end   = std::stod(currLineVec[11]);
  double qcov  = end - start;
  if (qcov < 0) qcov *= -1;
  qcov += 1;
  return qcov;
}

std::vector<std::string> expRange(std::vector<std::string> & hsp1,
                                  std::vector<std::string> & hsp2) {
  if (hsp1.empty()) {
    return hsp2;
  }
  if (hsp2.empty()) {
    return hsp1;
  }
  std::string start1 = hsp1[10];
  std::string end1   = hsp1[11];
  std::string start2 = hsp2[10];
  std::string end2   = hsp2[11];
  double start1d = std::stod(start1, NULL);
  double end1d   = std::stod(end1, NULL);
  double start2d = std::stod(start2, NULL);
  double end2d   = std::stod(end2, NULL);
  
  std::string start;
  std::string end;
  if (start1d < end1d && start2d < end2d) {
    start = (start1d < start2d) ? start1 : start2;
    end   = (end1d < end2d) ? end2 : end1;
  }
  else if (start1d > end1d && start2d > end2d) {
    start = (start1d > start2d) ? start1 : start2;
    end   = (end1d > end2d) ? end2 : end1;
  }
  else {
    return hsp1;
  }
  hsp1[10] = start;
  hsp1[11] = end;
  return hsp1;
}


bool isSeparate(std::vector<std::string> hsp1, std::vector<std::string> hsp2) {
  double len1 = getQcov(hsp1);
  double len2 = getQcov(hsp2);
  double start1 = std::stod(hsp1[10]);
  double end1   = std::stod(hsp1[11]);
  double start2 = std::stod(hsp2[10]);
  double end2   = std::stod(hsp2[11]);
  std::vector<double> bounds{start1, end1, start2, end2};
  double start = bounds[0];
  double end = bounds[0];
  for (auto mn : bounds) {
    if (mn < start) {
      start = mn;
    }
  }
  for (auto mx : bounds) {
    if (mx > end) {
      end = mx;
    }
  }
  double overlap = len1 + len2 - (end - start) + 1;
  double lenMin = (len1 < len2) ? len1 : len2;
  double ovMin = (60 < 0.2 * lenMin) ? 60 : 0.2 * lenMin;
  if (overlap < ovMin) {
    return true;
  }
  else {
    return false;
  }
}


bool checkBlock(std::vector<std::vector<std::string>> & block, bool multiGene,
                std::ofstream & outFile1, std::ofstream & outFile2) {
  if (block.size() == 1) {
    return true;
  }
  std::vector<std::string> pos;
  std::vector<std::string> neg;
  for (auto hsp : block) {
    if (hsp[4][0] == '-') {
      neg = expRange(neg, hsp);
    }
    else {
      pos = expRange(pos, hsp);
    }
  }
  if ((pos.empty() && !neg.empty()) || (!pos.empty() && neg.empty())) {
    return true;
  }
  else if (isSeparate(pos, neg)) {
    if (multiGene) {
      std::string posStart = pos[10];
      std::string posEnd   = pos[11];
      std::string negStart = neg[10];
      std::string negEnd   = neg[11];
      std::string start1 = (std::stod(posStart) < std::stod(posEnd)) ? posStart : posEnd;
      std::string end1   = (std::stod(posStart) > std::stod(posEnd)) ? posStart : posEnd;
      std::string start2 = (std::stod(negStart) < std::stod(negEnd)) ? negStart : negEnd;
      std::string end2   = (std::stod(negStart) > std::stod(negEnd)) ? negStart : negEnd;
      outFile1 << (pos[0] + " " + start1 + " " + end1 + " trans-multi\n");
      outFile1 << (neg[0] + " " + start2 + " " + end2 + " trans-multi\n");
    }
    else {
      std::vector<std::string> outHsp = (getQcov(pos) > getQcov(neg)) ? pos : neg;
      std::string start = (std::stod(outHsp[10]) < std::stod(outHsp[11])) ? outHsp[10] : outHsp[11];
      std::string end = (std::stod(outHsp[10]) > std::stod(outHsp[11])) ? outHsp[10] : outHsp[11];
      outFile1 << (outHsp[0] + " " + start + " " + end + " trans-self\n");
    }
    for (auto sp : pos) {
      outFile2 << sp + "\t";
    }
    outFile2 << "\n";
    for (auto sn : neg) {
      outFile2 << sn + "\t";
    }
    outFile2 << "\n";
  }
  else {
    return true;
  }
  return true;
}

void detect_chimera(std::string blastxFile, std::string outDir) {
  fs::path blastxFilePath(blastxFile.c_str());
  std::string blastxFileStr(blastxFilePath.stem().c_str());

  std::string cutFile = blastxFileStr + ".cut";
  std::string infoFile = blastxFileStr + ".info";

  fs::path cutFilePath((outDir + "/" + cutFile).c_str());
  fs::path infoFilePath((outDir + "/" + infoFile).c_str());
  if (fs::exists(cutFilePath) && fs::exists(infoFilePath)) {
    std::cout << "Chimera detection files found for: " + blastxFileStr << std::endl;
    return;
  }
  std::ifstream fileBlastx(blastxFile);
  std::ofstream fileCut(cutFilePath.c_str());
  std::ofstream fileInfo(infoFilePath.c_str());

  std::string currLine;
  std::vector<std::string> currLineVec;
  std::vector<std::vector<std::string>> hitB;
  std::vector<std::vector<std::string>> qryB;
  std::string query;
  std::string hit;
  std::string last_query;
  std::string last_hit;
  double pident;
  double qstart;
  double qend;
  double qcov;
  bool seqG;
  while (std::getline(fileBlastx, currLine)) {
    if (currLine.length() < 3) {
      continue;
    }
    currLineVec = getLineVec(currLine);
    pident = std::stod(currLineVec[5]);
    qcov = getQcov(currLineVec);
    
    if (pident < PIDENT_CUTOFF || qcov < LENGTH_CUTOFF) {
      continue;
    }
    query = currLineVec[0];
    hit   = currLineVec[2];
    if (last_query == "") {
      hitB.push_back(currLineVec);
      qryB.push_back(currLineVec);
      seqG = true;
    }
    else if (query == last_query) {
      qryB.push_back(currLineVec);
      if (seqG) {
        if (hit == last_hit) {
          hitB.push_back(currLineVec);
        }
        else {
          seqG = checkBlock(hitB, false, fileCut, fileInfo);
          hitB = {currLineVec};
        }
      }
    }
    else {
      if (seqG) {
        seqG = checkBlock(hitB, false, fileCut, fileInfo);
      }
      if (seqG) {
        seqG = checkBlock(qryB, true, fileCut, fileInfo);
      }
      qryB = {currLineVec};
      hitB = {currLineVec};
      seqG = true;
    }
    last_query = query;
    last_hit   = hit;
  }
  if (seqG) {
    seqG = checkBlock(hitB, false, fileCut, fileInfo);
  }
  if (seqG) {
    seqG = checkBlock(qryB, true, fileCut, fileInfo);
  }
  
  fileBlastx.close();
  fileCut.close();
  fileInfo.close();
}

// Create list of unique chimeric transcripts based on .cut file
std::set<std::string> makeChimeraSet(std::ifstream & chimFile) {
  std::set<std::string> chimSet;
  std::string currChim;
  std::string currLine;
  size_t wsPos = 0;
  while (std::getline(chimFile, currLine)) {
    wsPos = currLine.find(' ');
    if (wsPos != std::string::npos) {
      currChim = currLine.substr(0, wsPos);
      chimSet.insert(currChim);
    }
  }
  return chimSet;
}


void removeChimera(std::string transIn, std::string transOut,
                   std::string infoFilePath, std::string cutFilePath,
                   std::string ram_gb, std::string logFile) {
  fs::path transPath(transIn.c_str());
  std::string transPathStr(transPath.c_str());
  std::string transFileStr(transPath.stem().c_str());
  std::set<std::string> chimeraSet;

  std::string filtTrans = transOut;

  if (fs::exists(fs::path(filtTrans))) {
    logOutput("Filtered transcripts found for: " + transFileStr, logFile);
    return;
  }

  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;

  std::ifstream cutFile(cutFilePath);
  std::ifstream infoFile(infoFilePath);
  // Create unique list of chimeric transcript IDs
  chimeraSet = makeChimeraSet(cutFile);
  // Determine size of hash table
  uintmax_t numBytesTrans = fs::file_size(transPath);
  uintmax_t lenHashTable = numBytesTrans / 160;

  seqHash fastaHashTable(lenHashTable, transPath, ram_b);
  // Fill hash table with all transcripts
  // Iterate over chimera set
  //   For all sequences in chimera set
  //     Remove that sequence from the hash table
  //   Resulting hash table contains only non-chimeric seqs
  bool foundInHashTable;
  int numRemoved = 0;
  for (auto head : chimeraSet) {
    foundInHashTable = fastaHashTable.inHashTable(head);
    if (foundInHashTable) {
      fastaHashTable.deleteHash(head);
      foundInHashTable = fastaHashTable.inHashTable(head);
      if (foundInHashTable) {
        logOutput("Could not remove chimera. Something failed", logFile);
        return;
      }
      numRemoved++;
    }
    else {
      logOutput("Not found in hash table! Something failed", logFile);
      return;
    }
  }
  logOutput("Chimeric removal successful!\n", logFile);
  logOutput("====================", logFile);
  logOutput("  SUMMARY:", logFile);
  logOutput("    Removed:  " + std::to_string(numRemoved), logFile);
  logOutput("    Retained: " + std::to_string(fastaHashTable.getSize()), logFile);
  logOutput("====================\n", logFile);
  fastaHashTable.dump(filtTrans);
}

