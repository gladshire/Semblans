#include "annotate.h"

// Given a HMMER output file, produce a dictionary of query header - best hit
// header pairs
std::map<std::string, std::string> getGeneMatches(std::string annotFile) {

  std::map<std::string, std::string> geneMatches;
 
  if (fs::is_empty(fs::path(annotFile.c_str()))) {
    return geneMatches;
  }
  std::string transPrefix(fs::path(annotFile.c_str()).stem().stem().c_str());
  std::ifstream pantherFile(annotFile);
  std::string currLine;
  std::string currTok;
  std::string queryId;
  std::string lastQry;
  std::string matchId;
  std::string bestMatchId;
  std::string matchDesc;
  std::string bestMatchDesc;
  std::string eValStr;
  std::string lastEval;
  std::string newHeader;
  size_t tabPos;
  double eValCurr;
  double eValMin;
  int linNum = 0;

  // Parse tabular annotation file, capturing for each entry:
  //   - The header of the query
  //   - The PANTHER ID of the hit sequence
  //   - A functional description of the hit sequence
  //   - The hit's E-value
  //
  // Pair a query header with its best hit (the hit with the lowest E-value)
  while (getline(pantherFile, currLine)) {
    int colNum = 0;
    while ((tabPos = currLine.find("\t")) != std::string::npos) {
      currTok = currLine.substr(0, tabPos);
      if (colNum == 0) {
        queryId = currTok;
      }
      if (colNum == 1) {
        matchId = currTok;
      }
      if (colNum == 2) {
        matchDesc = currTok;
      }
      if (colNum == 3) {
        eValStr = currTok;
      }
      currLine.erase(0, tabPos + 1);
      colNum++;
    }
    std::istringstream doubleStr(eValStr);
    doubleStr >> eValCurr;
    if (queryId == lastQry) {
      if (eValCurr < eValMin) {
        doubleStr >> eValMin;
        bestMatchId = matchId;
        bestMatchDesc = matchDesc;
      }
    }
    else {
      doubleStr >> eValMin;
      bestMatchId = matchId;
      bestMatchDesc = matchDesc;
      if (lastQry != "") {
        newHeader = transPrefix + queryId.substr(queryId.find(" ")) + " " +
                    matchId + " " + matchDesc;
        geneMatches.emplace(queryId, newHeader);
      }
    }
    lastQry = queryId;
    lastEval = eValStr;
    linNum++;
  }
  return geneMatches;
}

// Perform annotation of transcripts using HMMER and PANTHER, generating a tabular output
// file, and then rename each transcript's header to one corresponding with its best hit
void annotateTranscript(std::string transIn, std::string transPep, std::string transOut,
                        std::string threads, std::string ram_gb, bool dispOutput,
                        std::string logFile) {
  fs::path transPepPath(transPep.c_str());
  std::string transPepFileStr(transPepPath.c_str());
  
  fs::path transInPath(transIn.c_str());
  std::string transInFileStr(transInPath.c_str());

  std::map<std::string, std::string> newSeqHeaders;
  // Prepare annotation data directory
  std::string annotFile((fs::path(transOut.c_str()).parent_path() / 
                         fs::path(transIn.c_str()).stem().stem()).c_str());
  annotFile += ".annot";

  // Determine size of hash table
  uintmax_t ram_b = (uintmax_t)stoi(ram_gb) * 1000000000;
  uintmax_t numBytesTrans = fs::file_size(transPepPath);
  uintmax_t lenHashTable = numBytesTrans / 160;

  // Create hash table with peptide sequences
  seqHash fastaPepHashTable(lenHashTable, transPepFileStr, ram_b);

  // Obtain data from hash table;
  //std::vector<sequence> * hashData = fastaPepHashTable.getHashData();
  linkedList * hashData = fastaPepHashTable.getHashData();

  // Initiate PANTHER scoring
  pantherScore(transPepFileStr, annotFile, threads, dispOutput, logFile);  

  // Obtain PANTHER gene descriptions
  newSeqHeaders = getGeneMatches(annotFile);

  // Create hash table with coding sequences
  seqHash fastaHashTable(lenHashTable, transInFileStr, ram_b);

  // Iterate over new seq ID data, renaming each header to new sequence ID
  for (auto seqPair : newSeqHeaders) {
    fastaHashTable.setSeqHeader(seqPair.first, seqPair.second);
  }

  // Dump annotated hash table to new location
  fastaHashTable.dump(transOut);
}
