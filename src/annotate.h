#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include "sra.h"
#include "transcript.h"
#include "seq.h"
#include "panther_score.h"
#include "seq_hash.h"
#include "print_info.h"

namespace fs = boost::filesystem;

std::map<std::string, std::string> getGeneMatches(std::string annotFile);

void annotateTranscript(std::string transIn, std::string transPep, std::string transOut,
                        std::string threads, std::string ram_gb, bool dispOutpu,
                        std::string logFile);
