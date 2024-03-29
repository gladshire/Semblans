#include "panther_score.h"


// Given a protein sequence, produce a scoring output file using HMMER against the
// PANTHER standard protein database
void pantherScore(std::string transPepIn, std::string outFile, std::string threads,
                  bool dispOutput, std::string logFile) {

  fs::path currDir = fs::current_path();
  fs::current_path(fs::path(PATH_PANTHER.c_str()).parent_path());
  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  if (fs::exists(outFile.c_str())) {
    return;
  }
  std::string panthCmd = "./pantherScore2.2.pl -V -l " + PANTHER_LIB + " -D A -H " +
                         PATH_HMMER + " -i " + transPepIn + " -o " + outFile + " -n -c " +
                         threads + printOut;

  result = system(panthCmd.c_str());
  fs::current_path(currDir);
  checkExitSignal(result, logFile);
}
