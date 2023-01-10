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
    std::string sequence;
    std::string quality;
    int numBp;
    int numLine;
  public:
    sequence();
    sequence(std::string header);
    sequence(std::string sequence);
    sequence(std::string header, std::string sequence);
    sequence(std::string sequence, std::string quality);
    sequence(std::string header, std::string sequence, std::string quality);
    seq
