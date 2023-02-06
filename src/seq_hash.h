#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "sra.h"
#include "transcript.h"
#include "seq.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;


class seqHash {
  private:
    std::vector<sequence> * seqHashData;
    uintmax_t lenHashTable;
  public:
    seqHash();
    seqHash(uintmax_t lenTable);
    unsigned long hashFunction(char * key);
    void insertHash(std::string header, std::string sequence);
    void deleteHash(std::string header);
    bool inHashTable(std::string header);
    void dump(std::string filePath);
};
