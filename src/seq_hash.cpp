#include "seq_hash.h"

// Default constructor for sequence hash table (seqHash) object
seqHash::seqHash() {
  seqHashData = nullptr;
  lenHashTable = 0;
  numItems = 0;
}

// Constructor for seqHash object allowing definition of hash table size
seqHash::seqHash(uintmax_t lenTable) {
  //seqHashData = new std::vector<sequence>[lenTable];
  seqHashData = new linkedList[lenTable]; 
  lenHashTable = lenTable;
  numItems = 0;
}

// Constructor for seqHash object providing FASTA/FASTQ file with sequences
// to fill it
seqHash::seqHash(uintmax_t lenTable, fs::path transFilePath, uintmax_t ram_b) {

  if (!fs::exists(transFilePath)) {
    seqHashData = nullptr;
    return;
  }
  //seqHashData = new std::vector<sequence>[lenTable];
  seqHashData = new linkedList[lenTable];
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
  char * qualityStartPos;
  char * qualityEndPos;
  char * inFileL;
  char headChar;
  while (!inFile.eof() && inFile.good()) {
    // Store file chunk into buffer
    inFile.read(&inFileData[0], ram_b);
    // Get number of bytes just read
    s = inFile.gcount();
    // Initialize header start position
    headerStartPos = &inFileData[0];
    if (*headerStartPos == '>') {
      headChar = '>';
    }
    else {
      headChar = '@';
    }
    // Initizlie address of last buffer position
    inFileL = &inFileData[0] + s;
    // Align end of buffer with end of last transcript
    this->align_buffer_end(inFile, &inFileData[0], s);
    std::string currHeader;
    std::string currSequence;
    std::string currQuality;
    std::string currKey;
    while (headerStartPos != inFileL) {
      // Extract header
      if (*headerStartPos != '>' && *headerStartPos != '@') {
        std::cout << "Buffer not aligned. Header doesn't start with \">\" or \"@\"" << std::endl;
        exit(0);
      }
      seqStartPos = std::find(headerStartPos, inFileL, '\n') + 1;
      currHeader = std::string(headerStartPos + 1, seqStartPos - 1);
      
      // Extract sequence
      seqEndPos = std::find(seqStartPos, inFileL, '\n');
      //currSequence = std::string(seqStartPos, seqEndPos);

      // Extract quality if FASTQ, insert sequence object into hash table
      if (*(seqEndPos + 1) == '+') {
        currSequence = std::string(seqStartPos, seqEndPos);
        qualityStartPos = std::find(seqEndPos + 1, inFileL, '\n') + 1;
        qualityEndPos = std::find(qualityStartPos, inFileL, '\n');
        if (qualityEndPos == inFileL) {
          currQuality = std::string(qualityStartPos, (int)(inFileL - qualityStartPos));
        }
        else {
          currQuality = std::string(qualityStartPos, qualityEndPos);
        }
        headerStartPos = std::find(qualityStartPos, inFileL, '\n') + 1;
        this->insertHash(currHeader, currSequence, currQuality);
      }
      else {
        headerStartPos = std::find(seqStartPos, inFileL, headChar);
        currSequence = std::string(seqStartPos, headerStartPos - 1);
        this->insertHash(currHeader, currSequence);
      }
    }
  }
  inFile.close();
}

// Utility function: remove characters from input stream until header character @ / > is reached
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
void seqHash::insertHash(std::string header, std::string sequenceData, std::string quality) {
  std::string keyStr = header.substr(0, header.find(' '));
  unsigned long hashIndex = hashFunction(&keyStr[0]) % lenHashTable;
  //seqHashData[hashIndex].push_back(sequence(header, sequenceData, quality));
  seqHashData[hashIndex].insert(sequence(header, sequenceData, quality));
  numItems++;
}

void seqHash::insertHash(std::string header, std::string sequenceData) {
  std::string keyStr = header.substr(0, header.find(' '));
  unsigned long hashIndex = hashFunction(&keyStr[0]) % lenHashTable;
  //seqHashData[hashIndex].push_back(sequence(header, sequenceData));
  seqHashData[hashIndex].insert(sequence(header, sequenceData));
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
  //auto dataIter = seqHashData[hashIndex].begin();
  /*
  while (dataIter != seqHashData[hashIndex].end()) {
    currKeyStrHash = (dataIter->get_header()).substr(0, dataIter->get_header().find(' '));
    if (keyStr == currKeyStrHash) {
      dataIter = seqHashData[hashIndex].erase(dataIter);
      numItems--;
    }
    else {
      dataIter++;
    }
  }*/
  bool removed = seqHashData[hashIndex].remove(header);
  if (removed) {
    numItems--;
  }
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
    //auto dataIter = seqHashData[hashIndex].begin();
    /*while (dataIter != seqHashData[hashIndex].end()) {
      currKeyStrHash = (dataIter->get_header()).substr(0, dataIter->get_header().find(' '));
      if (dataIter->get_header() == header || currKeyStrHash == keyStr) {
        return true;
      }
      dataIter++;
    }
    return false;
    */
    return seqHashData[hashIndex].exists(header);
  }
}

sequence seqHash::getSeq(std::string header) {
  std::string keyStr = header.substr(0, header.find(' '));
  unsigned long hashIndex = hashFunction(&keyStr[0]) % lenHashTable;
  if (seqHashData[hashIndex].empty()) {
    return sequence();
  }
  else {
    //std::string currKeyStrHash;
    //auto dataIter = seqHashData[hashIndex].begin();
    /*
    while (dataIter != seqHashData[hashIndex].end()) {
      currKeyStrHash = (dataIter->get_header()).substr(0, dataIter->get_header().find(' '));
      if (dataIter->get_header() == header || currKeyStrHash == keyStr) {
        return sequence(dataIter->get_header(), dataIter->get_sequence(),
                        dataIter->get_quality());
      }
      dataIter++;
    }
    return sequence();
    */
    //sequence seq = seqHashData[hashIndex].getSeq(header);
    return seqHashData[hashIndex].getSeq(header);
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
    //auto dataIter = seqHashData[hashIndex].begin();
    // Check for sequence in chained hash table
    /*
    while (dataIter != seqHashData[hashIndex].end()) {
      currKeyStrHash = (dataIter->get_header()).substr(0, dataIter->get_header().find(' '));
      // If match found update its header and return
      if (dataIter->get_header() == header || currKeyStrHash == keyStr) {
        dataIter->set_header(newHeader);
        return;
      }
      dataIter++;
    }
    */
    seqHashData[hashIndex].setSeqHead(header, newHeader);
    //std::cout << "ERROR: Sequence not found" << std::endl;
    //exit(2);
  }
}

// Output all hash table data to text file
void seqHash::dump(std::string filePath) {
  std::ofstream outFile(filePath);
  std::string currHeader;
  std::string currSeq;
  std::string currQual;
  for (uintmax_t i = 0; i < lenHashTable; i++) {
    if (seqHashData[i].empty()) {
      continue;
    }
    else {
      int bpPerLine;
      seqHashData[i].dump(outFile);
    }
  }
  outFile.close();
}

// Return length of hash table
uintmax_t seqHash::getLength() {
  return lenHashTable;
}

// Return number of items in hash table
uintmax_t seqHash::getNumItems() {
  return numItems;
}

// Get pointer to hash data
//std::vector<sequence> * seqHash::getHashData() {
linkedList * seqHash::getHashData() {
  return seqHashData;
}
/*
void seqHash::clear() {
  for (uintmax_t i = 0; i < lenHashTable; i++) {
    seqHashData[i].clear();
  }
  delete[] seqHashData;
}*/

// Destructor for hash table
seqHash::~seqHash() {
  if (seqHashData != nullptr) {
    for (uintmax_t i = 0; i < lenHashTable; i++) {
      seqHashData[i].clear();
    }
    delete[] seqHashData;
  }
}
