#pragma once
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>

typedef boost::iostreams::tee_device<std::ostream, std::ostream> teedev;
typedef boost::iostreams::stream<teedev, std::char_traits<typename std::ostream::char_type>,
                                 std::allocator<typename std::ostream::char_type>> teeStream;

void replaceChar(std::string inFilePath, char oldChar, char newChar);

void logOutput(std::string input, std::string logFile);

void printVertEllipse(std::string logFile, int numLines);

void printBreakLine(std::string logFile, int leadingSpace, int length);

std::string getPercent(float valPercent, int precision);

void reportError(int commandResult, std::string logFile);
