#include "ini_parse.h"


namespace fs = boost::filesystem;

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


// Create directory space for project
void make_proj_space(const INI_MAP &iniFile) {
  std::string projDir = iniFile.at("General").at("output_directory") +
                        iniFile.at("General").at("project_name") + "/";
  extern std::vector<std::string> stepDirs;
  system(("mkdir " + projDir).c_str());
  for (auto dir : stepDirs) {
    system(("mkdir " + projDir + dir).c_str());
  }
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

