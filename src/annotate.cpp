#include "annotate.h"


void annotateTranscript(std::transIn, std::string transOut,
                        std::string threads, std::string ram_gb,
                        std::string logFile, std::string email) {
  fs::path transPath(transIn.c_str());
  std::string transFileStr(transPath.stem().c_str());
  
  std::string annotatedTrans = transOut;

  // Check if checkpoint exists
  
  // Determine size of hash table
  uintmax_t numBytesTrans = fs::file_size(transPath);
  uintmax_t lenHashTable = numBytesTrans / 160;

  // Create hash table with sequences in transcript
  seqHash fastaHashTable(lenHashTable, transPath, ram_b);

  // Instantiate sequence ID job manager
  seqIdJobManager transAnnotator;

  // Obtain data from hash table;
  std::vector<sequence> * hashData = fastaHashTable.getHashData();

  // Iterate over hash table, submitting each sequence as job to annotator
  for (uintmax_t i = 0; i < lenHashTable; i++) {
    if (hashData[i].empty()) {
      continue;
    }
    else {
      for (auto seq : hashData[i]) {
        transAnnotator.submitSeqJob(seq);
      }
    }
  }

  // Initiate annotation manager
  transAnnotator.startSeqJobs(stoi(threads), email);

  // Obtain annotator new sequence ID data
  std::map<std::string, std::string> transAnnotator.getSeqIds();

  // Iterate over hash table, renaming each header to new sequence ID
  for (uintmax_t
}
