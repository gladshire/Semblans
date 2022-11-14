#include "rem_chimera.h"


void detect_chimera(std::string pathBlastxFile, std::string outDir) {
  fs::path blastxFilePath(pathBlastxFile.c_str());
  std::string blastxFileStr(blastxFilePath.stem().c_str());

  std::string cutFile = blastxFileStr + ".cut";
  std::string infoFile = blastxFileStr + ".info";

  fs::path cutFilePath((outDir + "/" + cutFile).c_str());
  fs::path infoFilePath((outDir + "/" + infoFile).c_str());
  if (fs::exists(cutFilePath) && fs::exists(infoFilePath)) {
    std::cout << "Chimera detection files found for: " + blastxFileStr;
    return;
  }

}

void rem_chimera(transcript transcripts, std::string infoFile);

void rem_chim_bulk(std::vector<transcript> transcriptsV, std::string infoFile);
