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
    std::string id;
    int numBp;
  public:
    sequence();
    sequence(std::string header, std::string sequenceData);
    sequence(const sequence & seq);
    std::string get_header();
    std::string get_sequence();
    std::string get_id();
    void set_id(std::string newId);
};
