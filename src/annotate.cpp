#include "annotate.h"


void annotateTranscript(std::string transIn, std::string transPep, std::string transOut,
                        std::string threads, std::string ram_gb,
                        std::string logFile, std::string email) {
  fs::path transPepPath(transPep.c_str());
  std::string transPepFileStr(transPepPath.c_str());
  
  fs::path transInPath(transIn.c_str());
  std::string transInFileStr(transInPath.c_str());

  std::string annotatedTrans = transOut;

  // Prepare annotation data directory
  std::string annotDir((fs::path(transOut.c_str()).parent_path() / 
                        fs::path(transIn.c_str()).stem().stem()).c_str());
  annotDir += ".annot_data";
  if (!fs::exists(fs::path(annotDir.c_str()))) {
    fs::create_directory(fs::path(annotDir.c_str()));
  }

  // Determine size of hash table
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  uintmax_t numBytesTrans = fs::file_size(transPepPath);
  uintmax_t lenHashTable = numBytesTrans / 160;

  // Create hash table with peptide sequences
  seqHash fastaPepHashTable(lenHashTable, transPepFileStr, ram_b);

  // Instantiate sequence ID job manager
  seqIdJobManager transAnnotator;

  // Obtain data from hash table;
  std::vector<sequence> * hashData = fastaPepHashTable.getHashData();

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
  transAnnotator.performSeqJobs(stoi(threads), email, annotDir);

  // Obtain annotator new sequence ID data
  std::cout << "Obtaining map with (oldSeq, newSeq) pairs" << std::endl;
  std::map<std::string, std::string> newSeqData = transAnnotator.getSeqIds();

  // Create hash table with coding sequences
  std::cout << "Creating coding seq hash table for transcript" << std::endl;
  seqHash fastaHashTable(lenHashTable, transInFileStr, ram_b);

  // Iterate over new seq ID data, renaming each header to new sequence ID
  std::cout << "Iterating over (oldSeq, newSeq) map, renaming headers to their new name" << std::endl;
  std::string currOldHeader;
  std::string currNewHeader;
  for (auto seqPair : newSeqData) {
    if (seqPair.second == "") {
      continue;
    }
    currOldHeader = seqPair.first;
    currNewHeader = std::string(fs::path(transOut.c_str()).stem().stem().c_str()) +
                    "_" + currOldHeader.substr(0, currOldHeader.find(" ")) + " " +
                    seqPair.second;
    fastaHashTable.setSeqHeader(seqPair.first, currNewHeader);
  }

  // Dump annotated hash table to new location
  std::cout << "Dumping annotated sequence hash to FASTA file" << std::endl;
  fastaHashTable.dump(transOut);
}
