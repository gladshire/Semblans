#include "seq_hash.h"


seqHash::seqHash() {
  seqHashData = NULL;
  lenHashTable = 0;
  numItems = 0;
}

seqHash::seqHash(uintmax_t lenTable) {
  seqHashData = new std::vector<sequence>[lenTable];
  lenHashTable = lenTable;
  numItems = 0;
}

seqHash::seqHash(uintmax_t lenTable, fs::path transFilePath, uintmax_t ram_b) {

  seqHashData = new std::vector<sequence>[lenTable];
  lenHashTable = lenTable;
  numItems = 0;

  std::string inFileStr(transFilePath.c_str());
  std::ifstream inFile(inFileStr);
  std::string inFileData;

  inFileData.reserve(ram_b);
  std::streamsize s;

  char * headerStartPos;
  char * headerEndPos;
  char * seqStartPos;
  char * seqEndPos;
  char * inFileL;
  while (!inFile.eof() && inFile.good()) {
    // Store file chunk into buffer
    inFile.read(&inFileData[0], ram_b);
    // Get number of bytes just read
    s = inFile.gcount();
    // Initialize header start position
    headerStartPos = &inFileData[0];
    // Initizlie address of last buffer position
    inFileL = &inFileData[0] + s;
    // Align end of buffer with end of last transcript
    this->align_buffer_end(inFile, &inFileData[0], s);
    std::string currHeader;
    std::string currSequence;
    std::string currKey;
    while (headerStartPos != inFileL) {
      // Extract header
      if (*headerStartPos != '>') {
        std::cout << "Buffer not aligned. Header doesn't start with \">\" or \"@\"" << std::endl;
        exit(0);
      }
      seqStartPos = std::find(headerStartPos, inFileL, '\n') + 1;
      currHeader = std::string(headerStartPos + 1, seqStartPos - 1);
      // Extract sequence
      headerStartPos = std::find(seqStartPos, inFileL, '>');
      currSequence = std::string(seqStartPos, headerStartPos - 1);
      // Insert information into hash table
      this->insertHash(currHeader, currSequence);
    }
  }
}

void seqHash::align_buffer_end(std::ifstream & inFile, char * inFileData, std::streamsize & s) {
  if (!inFile.eof() && inFile.good()) {
    while (inFile.peek() != '@' && inFile.peek() != '>') {
      std::cout << inFile.peek() << std::endl;
      inFile.unget();
      inFileData[s - 1] = '\0';
      s--;
    }
  }
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
  numItems++;
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
  numItems--;
}

// Determine if sequence element is contained in hash table
bool seqHash::inHashTable(std::string header) {
  std::string keyStr = header.substr(0, header.find(' '));
  unsigned long hashIndex = hashFunction(&keyStr[0]) % lenHashTable;
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

sequence seqHash::getSeq(std::string header) {
  std::string keyStr = header.substr(0, header.find(' '));
  unsigned long hashIndex = hashFunction(&keyStr[0]) % lenHashTable;
  if (seqHashData[hashIndex].empty()) {
    return sequence();
  }
  else {
    std::string currKeyStrHash;
    auto vecIter = seqHashData[hashIndex].begin();
    while (vecIter != seqHashData[hashIndex].end()) {
      currKeyStrHash = (vecIter->get_header()).substr(0, vecIter->get_header().find(' '));
      if (vecIter->get_header() == header || currKeyStrHash == keyStr) {
        return sequence(vecIter->get_header(), vecIter->get_sequence());
      }
      vecIter++;
    }
    return sequence();
  }
}

// Set new sequence header
void seqHash::setSeqHeader(std::string header, std::string newHeader) {
  std::string keyStr = header.substr(0, header.find(' '));
  unsigned long hashIndex = hashFunction(&keyStr[0]) % lenHashTable;
  if (seqHashData[hashIndex].empty()) {
    std::cout << "ERROR: Sequence not found" << std::endl;
    exit(2);
  }
  else {
    std::string currKeyStrHash;
    auto vecIter = seqHashData[hashIndex].begin();
    // Check for sequence in chained hash table
    while (vecIter != seqHashData[hashIndex].end()) {
      // If match found update its header and return
      if (vecIter->get_header() == header || currKeyStrHash == keyStr) {
        vecIter->set_header(newHeader);
        return;
      }
    }
    std::cout << "ERROR: Sequence not found" << std::endl;
    exit(2);
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
      int bpPerLine;
      for (auto seq : seqHashData[i]) {
        currHeader = seq.get_header();
        currSeq = seq.get_sequence();
        outFile << ">" << currHeader << '\n';
        outFile << currSeq << std::endl;
      }
    }
  }
  outFile.close();
}

// Return number of items in hash table
uintmax_t seqHash::getSize() {
  return numItems;
}

// Get pointer to hash data
std::vector<sequence> * seqHash::getHashData() {
  return seqHashData;
}
