#include "annotate.h"


void annotateTranscript(std::string transIn, std::string transPep, std::string transOut,
                        std::string threads, std::string ram_gb,
                        std::string logFile, std::string email) {
  fs::path transPepPath(transPep.c_str());
  std::string transPepFileStr(transPepPath.c_str());
  
  fs::path transInPath(transIn.c_str());
  std::string transInFileStr(transInPath.c_str());

  std::string annotatedTrans = transOut;

  // Determine size of hash table
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  uintmax_t numBytesTrans = fs::file_size(transPepPath);
  uintmax_t lenHashTable = numBytesTrans / 160;

  // Create hash table with peptide sequences
  std::cout << "Creating peptide hash table for transcript" << std::endl;
  seqHash fastaPepHashTable(lenHashTable, transPepFileStr, ram_b);

  // Instantiate sequence ID job manager
  std::cout << "Instantiating job manager for annotation" << std::endl;
  seqIdJobManager transAnnotator;

  // Obtain data from hash table;
  std::cout << "Obtaining peptide sequence data from hash table" << std::endl;
  std::vector<sequence> * hashData = fastaPepHashTable.getHashData();

  // Iterate over hash table, submitting each sequence as job to annotator
  std::cout << "Iterating over hash table, submitting annotation jobs for seqs" << std::endl;
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
  std::cout << "Starting sequence annotation jobs" << std::endl;
  transAnnotator.startSeqJobs(stoi(threads), email);

  // Obtain annotator new sequence ID data
  std::cout << "Obtaining map with (oldSeq, newSeq) pairs" << std::endl;
  std::map<std::string, std::string> newSeqData = transAnnotator.getSeqIds();

  // Create hash table with coding sequences
  std::cout << "Creating coding seq hash table for transcript" << std::endl;
  seqHash fastaHashTable(lenHashTable, transInFileStr, ram_b);

  // Iterate over new seq ID data, renaming each header to new sequence ID
  std::cout << "Iterating over (oldSeq, newSeq) map, renaming headers to their new name" << std::endl;
  for (auto seqPair : newSeqData) {
    fastaHashTable.setSeqHeader(seqPair.first, seqPair.second);
  }

  // Dump annotated hash table to new location
  std::cout << "Dumping annotated sequence hash to FASTA file" << std::endl;
  fastaHashTable.dump(transOut);
}
