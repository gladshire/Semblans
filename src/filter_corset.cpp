#include "filter_corset.h"

void fillFastaHash(seqHash fastaHashTable, transcript trans, uintmax_t ram_b) {
  fs::path transFilePathC = trans.get_trans_path_chim_filt();

  std::string inFileStr(transFilePathC.c_str());

  std::ifstream inFile(inFileStr);

  std::string inFileData;

  inFileData.reserve(ram_b);

  std::streamsize s;

  char * headerStartPos;
  char * headerEndPos;
  char * seqStartPos;
  char * seqEndPos;
  char * inFileL;
  while (!inFile.eof()) {
    // Process file in buffer chunks
    // Store chunk into buffer
    inFile.read(&inFileData[0], ram_b);
    // Get number of bytes just read
    s = inFile.gcount();
    // Initialize header start position
    headerStartPos = &inFileData[0];
    // Initizlie address of last buffer position
    inFileL = &inFileData[0] + s;
    // Align end of buffer just before next transcript
    align_buffer_end(inFile, &inFileData[0], s);
    // While loop to retrieve all transcripts in buffer
    //   Each iteration corresponds to single transcript
    std::string currHeader;
    std::string currSequence;
    std::string currKey;
    while (headerStartPos != inFileL) {
      // Extract header
      if (*headerStartPos != '>') {
        std::cout << "Buffer not aligned. Header doesn't start with \">\"" << std::endl;
        exit(0);
      }
      seqStartPos = std::find(headerStartPos, inFileL, '\n') + 1;
      currHeader = std::string(headerStartPos + 1, seqStartPos - 1);
      // Extract sequence
      headerStartPos = std::find(seqStartPos, inFileL, '>');
      currSequence = std::string(seqStartPos, headerStartPos - 1);
      // Insert information into hash table
      fastaHashTable.insertHash(currHeader, currSequence);
    }
  }
}

std::set<std::string> makeClusterSet(std::ifstream & clustFile) {
  std::set<std::string> clustSet;
  std::string currClust;
  std::string currLine;
  size_t wsPos = 0;
  while (std::getline(clustFile, currLine)) {
    // Find position of tab or space
    // wsPos = ??
    if (wsPos != std::string::npos) {
      currClust = currLine.substr(0, wsPos);
      clustSet.insert(currClust);
    }
  }
  return clustSet;
}

void filterCorset(transcript trans, std::string clusterPath, uintmax_t ram_b,
                  std::string out_dir) {
  fs::path transPath = trans.get_trans_path_chimera();
  std::string clustPath(trans.get_trans_path_clust().c_str());
  std::string largestClust(trans.get_trans_path_largest().c_str());
  std::string redundTrans(trans.get_trans_path_redund().c_str());

  if (fs::exists(largestClust) && fs::exists(redundTrans)) {
    std::cout << "Largest cluster and redundant transcripts found for: " << transPath
              << std::endl;
    return;
  }
  // Read largest corset clusters into set
  // Read transcripts into hash table
  // For each item in cluster set
  //   Retrieve seq from hash table
  //   Place into new hash table
  //   Remove from original hash table
  std::ifstream clustFile(clusterPath);

  uintmax_t numBytesTrans = fs::file_size(transPath);
  uintmax_t lenHashTable = numBytesTrans / 160;
  uintmax_t hashIndex;
  std::set<std::string> clustSet;

  clustSet = makeClusterSet(clustFile);
  seqHash seqRemovedHash(lenHashTable);
  seqHash seqRetainedHash(lenHashTable);

  fillFastaHash(seqRemovedHash, trans, ram_b);
  bool foundInHashTable;
  int numRemoved = 0;
  sequence currSeq;
  for (auto head : clustSet) {
    foundInHashTable = seqRemovedHash.inHashTable(head);
    if (foundInHashTable) {
      // Insert seq into seqRetainedHash
      currSeq = seqRemovedHash.getSeq(head);
      seqRetainedHash.insertHash(currSeq.get_header(), currSeq.get_sequence());
      // Delete seq from seqRemovedHash
      seqRemovedHash.deleteHash(head);
    }
  }
  // Dump largest cluster data to new fasta file
  seqRetainedHash.dump(largestClust); 
}
