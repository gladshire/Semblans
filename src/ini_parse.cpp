#include "ini_parse.h"


// Define output directory names in steps of pipeline
std::vector<std::string> allSteps = {"Raw_reads/", "Pre_quality_analysis/",
                                     "Error_correction/", "Trimming/",
                                     "Filter_foreign/", "Quality_analysis/",
                                     "Filter_overrepresented/", "Transcript_assembly/",
                                     "Filter_chimera/", "Clustering/",
                                     "Coding_seqs_peps/", "Final_annotated/"};

std::vector<std::string> stepDirs;


namespace fs = boost::filesystem;


// Convert stringified INI array to string vector
std::vector<std::string> getStrArray(std::string iniArrStr, std::string delim) {
  std::vector<std::string> strArray;
  if (iniArrStr.empty()) {
    return strArray;
  }
  size_t delimPos = 0;
  std::string currSubstr;
  while ((delimPos = iniArrStr.find(delim)) != std::string::npos) {
    currSubstr = iniArrStr.substr(0, delimPos);
    strArray.push_back(currSubstr);
    iniArrStr.erase(0, delimPos + delim.length());
  }
  strArray.push_back(iniArrStr);
  return strArray;
}

// Dispatch INI data to a map, indexable by sections, then by keys
static int ini_callback(IniDispatch * const dispatch, void * map_pt) {
  #define thismap (reinterpret_cast<INI_MAP*>(map_pt))
  if (dispatch->type == INI_COMMENT) {
    return 0;
  }
  if (dispatch->type == INI_SECTION) {
    INI_MAP_ENTRY newSec;
    thismap->insert(std::pair<std::string, INI_MAP_ENTRY>(dispatch->data, newSec));
    return 0;
  }
  if (dispatch->type == INI_KEY) {
    (*thismap)[dispatch->append_to].insert(std::pair<std::string, std::string>(dispatch->data, dispatch->value));
  }
  return 0;
}

std::string getStepNum(int stepNum) {
  std::string currNum;
  currNum = std::to_string(stepNum);
  if (stepNum < 10) {
    currNum = "0" + currNum;
  }
  return currNum;
}

// Create directory space for project
void make_proj_space(const INI_MAP &iniFile, std::string pipeStage) {
  std::string projDir((fs::path(iniFile.at("General").at("output_directory").c_str()) /
                       fs::path(iniFile.at("General").at("project_name").c_str())).c_str());
  projDir += "/";
  INI_MAP_ENTRY pipeSteps = iniFile.at("Pipeline");
  if (!fs::exists(fs::path(projDir.c_str()) / ".checkpoints")) {
    system(("mkdir " + projDir + "/.checkpoints/").c_str());
  }

  int stepNum = 0;
  std::string currNum = "00";

  // Create directory for raw reads
  stepDirs.push_back(currNum + "-" + allSteps[0]);
  system(("mkdir " + projDir + currNum + "-" + allSteps[0] + " > /dev/null 2>&1").c_str());
  stepNum++;

  // Create directory for initial quality analysis by FastQC
  currNum = getStepNum(stepNum);
  stepDirs.push_back(currNum + "-" + allSteps[1]);
  system(("mkdir " + projDir + currNum + "-" + allSteps[1] + " > /dev/null 2>&1").c_str());
  stepNum++;
  
  // If pipeline running error correction, create its directory
  currNum = getStepNum(stepNum);
  if (ini_get_bool(pipeSteps.at("error_correction").c_str(), -1)) {
    stepDirs.push_back(currNum + "-" + allSteps[2]);
    system(("mkdir " + projDir + currNum + "-" + allSteps[2] + " > /dev/null 2>&1").c_str());
    stepNum++;
  }
  else {
    stepDirs.push_back(allSteps[2]);
  }

  // If pipeline running adapter trimming, create its directory
  currNum = getStepNum(stepNum);
  if (ini_get_bool(pipeSteps.at("trim_adapter_seqs").c_str(), -1)) {
    stepDirs.push_back(currNum + "-" + allSteps[3]);
    system(("mkdir " + projDir + currNum + "-" + allSteps[3] + " > /dev/null 2>&1").c_str());
    stepNum++;
  }
  else {
    stepDirs.push_back(allSteps[3]);
  }

  // If pipeline running foreign sequence filtering, create its directory
  currNum = getStepNum(stepNum);
  if (ini_get_bool(pipeSteps.at("filter_foreign_reads").c_str(), -1)) {
    stepDirs.push_back(currNum + "-" + allSteps[4]);
    system(("mkdir " + projDir + currNum + "-" + allSteps[4] + " > /dev/null 2>&1").c_str());
    stepNum++;
  }
  else {
    stepDirs.push_back(allSteps[4]);
  }

  // Create directory for second quality check
  currNum = getStepNum(stepNum);
  stepDirs.push_back(currNum + "-" + allSteps[5]);
  system(("mkdir " + projDir + currNum + "-" + allSteps[5] + " > /dev/null 2>&1").c_str());
  stepNum++;

  // If pipeline running removal of overrepresented, create its directory
  currNum = getStepNum(stepNum);
  if (ini_get_bool(pipeSteps.at("remove_overrepresented").c_str(), -1)) {
    stepDirs.push_back(currNum + "-" + allSteps[6]);
    system(("mkdir " + projDir + currNum + "-" + allSteps[6] + " > /dev/null 2>&1").c_str());
    stepNum++;
  }
  else {
    stepDirs.push_back(allSteps[6]);
  }

  // Create directory for transcriptome assembly
  currNum = getStepNum(stepNum);
  stepDirs.push_back(currNum + "-" + allSteps[7]);
  system(("mkdir " + projDir + currNum + "-" + allSteps[7] + " > /dev/null 2>&1").c_str());
  stepNum++;

  // If pipeline running removal of chimeras, create its directory
  currNum = getStepNum(stepNum);
  if (ini_get_bool(pipeSteps.at("remove_chimera_reads").c_str(), -1)) {
    stepDirs.push_back(currNum + "-" + allSteps[8]);
    system(("mkdir " + projDir + currNum + "-" + allSteps[8] + " > /dev/null 2>&1").c_str());
    stepNum++;
  }
  else {
    stepDirs.push_back(allSteps[8]);
  }
  
  // If pipeline running cluster-based filtration, create its directory
  currNum = getStepNum(stepNum);
  if (ini_get_bool(pipeSteps.at("cluster_filtering").c_str(), -1)) {
    stepDirs.push_back(currNum + "-" + allSteps[9]);
    system(("mkdir " + projDir + currNum + "-" + allSteps[9] + " > /dev/null 2>&1").c_str());
    stepNum++;
    }
  else {
    stepDirs.push_back(allSteps[9]);
  }
  
  // Create directory for transcripts coding seqs, proteins
  currNum = getStepNum(stepNum);
  stepDirs.push_back(currNum + "-" + allSteps[10]);
  system(("mkdir " + projDir + currNum + "-" + allSteps[10] + " > /dev/null 2>&1").c_str());
  stepNum++;

  // Create directory for annotated transcripts
  currNum = getStepNum(stepNum);
  stepDirs.push_back(currNum + "-" + allSteps[11]);
  system(("mkdir " + projDir + currNum + "-" + allSteps[11] + " > /dev/null 2>&1").c_str());
  stepNum++;
}


// Given the path/name of an INI config file, return a map of its data
// which can be indexed by sections, and then by keys
INI_MAP make_ini_map(const char * configPath) {
  FILE * configIni = fopen(configPath, "r");
  fs::path configPathObj(configPath);
  try {
    if (!configIni) {
      std::string fileError = "ERROR: Cannot open config file: ";
      throw std::runtime_error(fileError);
    }
  } catch (std::runtime_error& e){
    std::cerr << e.what() << configPathObj.filename() << std::endl;
    return {};
  }
  INI_MAP iniMap;
  load_ini_file(configIni, INI_DEFAULT_FORMAT, NULL, ini_callback, &iniMap);
  fclose(configIni);
  return iniMap;
}

