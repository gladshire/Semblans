#pragma once
#include <fstream>
#include <iostream>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include "log.h"
#include "sra.h"
#include "transcript.h"

extern std::atomic<bool> procRunning;

void progressAnim(std::string precedeString, std::string logFile);

void summarize_sing_sra(SRA sra, std::string logFile,
                        int margin);

void summarize_all_sras(const std::vector<SRA> & sras, std::string logFile,
                        int margin);

void summarize_sing_trans(transcript trans, std::string logFile, int margin);

void summarize_all_trans(const std::vector<transcript> & transVec,
                         std::string logFile, int margin);
