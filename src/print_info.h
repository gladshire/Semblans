#pragma once
#include <fstream>
#include <iostream>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include "sra.h"
#include "transcript.h"

typedef boost::iostreams::tee_device<std::ostream, std::ostream> teedev;
typedef boost::iostreams::stream<teedev, std::char_traits<typename std::ostream::char_type>,
                                 std::allocator<typename std::ostream::char_type>> teeStream;

void logOutput(std::string input, std::string logFile);

void summarize_sing_sra(SRA sra, std::string logFile,
                        int margin);

void summarize_all_sras(const std::vector<SRA> & sras, std::string logFile,
                        int margin);

void summarize_sing_trans(transcript trans, std::string logFile, int margin);

void summarize_all_trans(const std::vector<transcript> & transVec,
                         std::string logFile, int margin);
