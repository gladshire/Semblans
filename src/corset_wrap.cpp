#include "corset_wrap.h"

void corset_eq_classes(std::string transPrefix, std::string pathEqClassFile,
                       std::string outDir, bool dispOutput, std::string logFile) {
  fs::path eqClassFile(pathEqClassFile);
  if (eqClassFile.extension() == fs::path(".gz")) {
    system(("gzip -d " + pathEqClassFile).c_str());
    pathEqClassFile = std::string(eqClassFile.replace_extension().c_str());
  }
  std::string printOut;
  if (dispOutput) {
    printOut = " 2>&1 | tee -a " + logFile;
  }
  else {
    printOut = " >>" + logFile + " 2>&1";
  }
  int result;
  std::string cors_cmd = PATH_CORSET + " -i" + " salmon_eq_classes " +
                         pathEqClassFile + " -m" + " 5" +
                         " -p " + outDir + "/" + transPrefix + "_salmon" +
                         printOut;
  result = system(cors_cmd.c_str());
  if (WIFSIGNALED(result)) {
    logOutput("Exited with signal " + std::to_string(WTERMSIG(result)), logFile);
    exit(1);
  }
}
