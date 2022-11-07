// Header file for ini_parse.cpp
// Enables the parsing of INI config files into the data type looking like this:
//   std::map<std::string, std::map<std::string, std::string>>



#include <iostream>
#include <string>
#include <map>
#include <cstdio>

#include <confini.h>
#include <boost/filesystem.hpp>
#include <boost/variant.hpp>

// Create short aliases which represent the map for the INI file
#define INI_MAP_ENTRY std::map<std::string, std::string>
#define INI_MAP std::map<std::string, INI_MAP_ENTRY>



static int ini_callback(IniDispatch * const dispatch, void * map_pt);

void make_proj_space(const INI_MAP &iniFile);

INI_MAP make_ini_map(const char * configPath);

