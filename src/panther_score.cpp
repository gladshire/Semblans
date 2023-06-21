#include "panther_score.h"

void pantherScore(std::string transPepIn, std::string outFile, std::string threads) {

  std::string panthCmd = PATH_PANTHER + " -l " + PANTHER_LIB + " -D B -H " + PATH_HMMER +
                         " -i " + transPepIn + " -o " + outFile + " -n -c " + threads;
  system(panthCmd.c_str());

}
