#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <functional>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "sra.h"
#include "transcript.h"
#include "seq.h"
#include "llist.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;


class seqHash {
  private:
    //std::vector<sequence> * seqHashData;
    linkedList * seqHashData;
    uintmax_t lenHashTable;
    uintmax_t numItems;
  public:
    seqHash();
    seqHash(uintmax_t lenTable);
    seqHash(uintmax_t lenTable, fs::path transFilePath, uintmax_t ram_b);
    void align_buffer_end(std::ifstream & inFile, char * inFileData, std::streamsize & s);
    unsigned long hashFunction(char * key);
    void insertHash(std::string header, std::string sequence, std::string quality);
    void insertHash(std::string header, std::string sequence);
    void deleteHash(std::string header);
    bool inHashTable(std::string header);
    sequence getSeq(std::string header);
    void setSeqHeader(std::string header, std::string newHeader);
    void dump(std::string filePath);
    uintmax_t getLength();
    uintmax_t getNumItems();
    void clear();
    linkedList * getHashData();    
    ~seqHash();
};
