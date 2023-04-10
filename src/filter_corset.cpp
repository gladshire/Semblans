#include "filter_corset.h"


std::set<std::string> makeClusterSet(std::ifstream & clustFile) {
  std::set<std::string> clustSet;
  std::string currClust;
  std::string currLine;
  size_t wsPos = 0;
  while (std::getline(clustFile, currLine)) {
    // Find position of tab or space
    wsPos = currLine.find_first_of("\t ");
    if (wsPos != std::string::npos) {
      currClust = currLine.substr(0, wsPos);
      clustSet.insert(currClust);
    }
  }
  return clustSet;
}

void filterCorset(std::string transIn, std::string transClust,
                  std::string transLargestClust, std::string transRedund,
                  uintmax_t ram_b, std::string outDir, std::string logFile) {
  fs::path transPath(transIn.c_str());
  std::string largestClust(transLargestClust);
  std::string redundTrans(transRedund);

  if (fs::exists(fs::path(largestClust)) && fs::exists(fs::path(redundTrans))) {
    logOutput("Largest cluster and redundant transcripts found for: " +
              std::string(transPath.c_str()), logFile);
    return;
  }
  // Read largest corset clusters into set
  // Read transcripts into hash table
  // For each item in cluster set
  //   Retrieve seq from hash table
  //   Place into new hash table
  //   Remove from original hash table
  std::ifstream clustFile(transClust);

  uintmax_t numBytesTrans = fs::file_size(transPath);
  uintmax_t lenHashTable = numBytesTrans / 160;
  uintmax_t hashIndex;
  std::set<std::string> clustSet;

  clustSet = makeClusterSet(clustFile);
  seqHash seqRemovedHash(lenHashTable, transPath, ram_b);
  seqHash seqRetainedHash(lenHashTable);

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
    else {
      logOutput("Not found in hash table! Something failed", logFile);
      return;
    }
  }
  // Dump largest cluster data to new fasta file
  seqRetainedHash.dump(largestClust); 
}
