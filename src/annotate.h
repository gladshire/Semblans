#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include "sra.h"
#include "transcript.h"
#include "seq.h"
#include "ips_job_man.h"
#include "seq_hash.h"


namespace fs = boost::filesystem;


void annotateTranscript(std::string transIn, std::string transPep, std::string transOut,
                        std::string threads, std::string ram_gb,
                        std::string logFile, std::string email);
