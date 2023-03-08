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

void filterCorset(transcript trans, std::string clusterPath, uintmax_t ram_b,
                  std::string out_dir, std::string logFile) {
  fs::path transPath = trans.get_trans_path_chimera();
  std::string clustPath(trans.get_trans_path_clust().c_str());
  std::string largestClust(trans.get_trans_path_largest().c_str());
  std::string redundTrans(trans.get_trans_path_redund().c_str());

  if (fs::exists(largestClust) && fs::exists(redundTrans)) {
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
  std::ifstream clustFile(clusterPath);

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
