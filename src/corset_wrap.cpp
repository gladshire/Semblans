#include "corset_wrap.h"

void corset_eq_classes(transcript trans, std::string pathEqClassFile, std::string outDir,
                       bool dispOutput, std::string logFile) {
  fs::path eqClassFile(pathEqClassFile);
  if (fs::exists(trans.get_trans_path_clust()) &&
      fs::exists(trans.get_trans_path_counts())) {
    std::cout << "Cluster/count files found for: " << trans.get_org_name() << std::endl;
    return;
  }
  if (eqClassFile.extension() == fs::path(".gz")) {
    system(("gzip -d " + pathEqClassFile).c_str());
    pathEqClassFile = std::string(eqClassFile.replace_extension().c_str());
  }
  std::string printOut;
  if (dispOutput) {
    printOut = " |& tee -a " + logFile;
  }
  else {
    printOut = " &>> " + logFile;
  }
  int result;
  std::string cors_cmd = PATH_CORSET + " -i" + " salmon_eq_classes " +
                         pathEqClassFile + " -m" + " 5" +
                         " -p " + outDir + trans.make_file_str() + "_salmon" +
                         printOut;
  result = system(cors_cmd.c_str());
  if (WIFSIGNALED(result)) {
    std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
    exit(1);
  }
}
