#pragma once
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>

typedef boost::iostreams::tee_device<std::ostream, std::ostream> teedev;
typedef boost::iostreams::stream<teedev, std::char_traits<typename std::ostream::char_type>,
                                 std::allocator<typename std::ostream::char_type>> teeStream;

void replaceChar(std::string inFilePath, char oldChar, char newChar);

void logOutput(std::string input, std::string logFile);

std::string getPercent(float valPercent, int precision);
