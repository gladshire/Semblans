#include "seq_hash.h"


seqHash::seqHash() {
  seqHashData = NULL;
  lenHashTable = 0;
}

seqHash::seqHash(uintmax_t lenTable) {
  seqHashData = new std::vector<sequence>[lenTable];
  lenHashTable = lenTable;
}

// Apply djb2 hashing function to key
unsigned long seqHash::hashFunction(char * key) {
  unsigned long hash = 5381;
  int c;
  while (c = *key++) {
    hash = ((hash << 5) + hash) + c;
  }
  return hash;
}

// Insert sequence element into hash table
void seqHash::insertHash(std::string header, std::string sequenceData) {
  std::string keyStr = header.substr(0, header.find(' '));
  unsigned long hashIndex = hashFunction(&keyStr[0]) % lenHashTable;
  seqHashData[hashIndex].push_back(sequence(header, sequenceData));
}

// Remove sequence element from hash table, if it exists
void seqHash::deleteHash(std::string header) {
  std::string keyStr = header.substr(0, header.find(' '));
  unsigned long hashIndex = hashFunction(&keyStr[0]) % lenHashTable;
  if (seqHashData[hashIndex].empty()) {
    return;
  }
  std::string currKeyStrHash;
  auto vecIter = seqHashData[hashIndex].begin();
  while (vecIter != seqHashData[hashIndex].end()) {
    currKeyStrHash = (vecIter->get_header()).substr(0, vecIter->get_header().find(' '));
    if (keyStr == currKeyStrHash) {
      vecIter = seqHashData[hashIndex].erase(vecIter);
    }
    else {
      vecIter++;
    }
  }
}

// Determine if sequence element is contained in hash table
bool seqHash::inHashTable(std::string header) {
  std::string keyStr = header.substr(0, header.find(' '));
  unsigned long hashIndex = hashFunction(&keyStr[0]) % lenHashTable;
  std::cout << "Index: " << hashIndex << std::endl;
  if (seqHashData[hashIndex].empty()) {
    return false;
  }
  else {
    std::string currKeyStrHash;
    auto vecIter = seqHashData[hashIndex].begin();
    while (vecIter != seqHashData[hashIndex].end()) {
      currKeyStrHash = (vecIter->get_header()).substr(0, vecIter->get_header().find(' '));
      if (vecIter->get_header() == header || currKeyStrHash == keyStr) {
        return true;
      }
      vecIter++;
    }
    return false;
  }
}

// Output all hash table data to text file
void seqHash::dump(std::string filePath) {
  std::ofstream outFile(filePath);
  std::string currHeader;
  std::string currSeq;
  for (uintmax_t i = 0; i < lenHashTable; i++) {
    if (seqHashData[i].empty()) {
      continue;
    }
    else {
      for (auto seq : seqHashData[i]) {
        currHeader = seq.get_header();
        currSeq = seq.get_sequence();
        outFile << currHeader << '\n';
        outFile << currSeq << '\n';
      }
    }
  }
  outFile.close();
}
