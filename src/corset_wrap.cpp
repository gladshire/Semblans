#include "corset_wrap.h"

// Given the output equivalence class file from a Salmon quant run, cluster transcript
// contigs into genes and identify redundant transcripts using Corset
void corset_eq_classes(std::string transPrefix, std::string pathEqClassFile,
                       std::string outDir, bool dispOutput, std::string logFile) {
  fs::path eqClassFile(pathEqClassFile);
  if (eqClassFile.extension() == fs::path(".gz")) {
    system(("gzip -d " + pathEqClassFile).c_str());
    pathEqClassFile = std::string(eqClassFile.replace_extension().c_str());
  }
  int result;
  std::string corsCmd = PATH_CORSET + " -i" + " salmon_eq_classes " +
                        pathEqClassFile + " -m" + " 5" +
                        " -p " + outDir + "/" + transPrefix + "_salmon";
  if (dispOutput) {
    corsCmd += " 2>&1 | tee -a " + logFile;
    logOutput("  Running command: " + corsCmd + "\n\n", logFile);
  }
  else {
    corsCmd += " >>" + logFile + " 2>&1";
  }
  result = system(corsCmd.c_str());
  checkExitSignal(result, logFile);
}
