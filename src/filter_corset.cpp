#include "filter_corset.h"

void filterCorset(transcript trans, std::string clusterPath, std::string out_dir) {
  std::string transPath(trans.get_trans_path_chimera().c_str());
  std::string largestClust(trans.get_trans_path_largest().c_str());
  std::string redundTrans(trans.get_trans_path_redund().c_str());

  if (fs::exists(largestClust) && fs::exists(redundTrans)) {
    std::cout << "Largest cluster and redundant transcripts found for: " << transPath
              << std::endl;
    return;
  }
  // Read corset clusters into dataframe

}
