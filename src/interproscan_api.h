#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <chrono>
#include <thread>
#include <boost/filesystem.hpp>
#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_utils.hpp>
#include <curl/curl.h>
#include "ini_parse.h"
#include "sra.h"
#include "transcript.h"
#include "seq.h"
#include "seq_hash.h"
#include "print_info.h"


// Define URL for interproscan API
#define IPS_URL std::string("https://www.ebi.ac.uk/Tools/services/rest/iprscan5")


size_t curlWriter(void * buffer, size_t size, size_t nmemb, std::string * s);

std::string submitJob(std::string email, std::string title, std::string sequence);

//std::string statusJob(std::string jobID);

//std::string resultType(std::string jobID);


