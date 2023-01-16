#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <thread>
#include <boost/filesystem.hpp>
#include "ini_parse.h"

namespace fs = boost::filesystem;

class sequence {
  private:
    std::string header;
    std::string sequenceData;
    std::string quality;
    int numBp;
    int numLine;
    int bpPerLine;
  public:
    sequence();
    sequence(std::string header, std::string sequenceData);
    std::string get_header();
    std::string get_sequence();
};
