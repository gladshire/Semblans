#include <string>
#include <vector>
#include <iostream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <curl/curl.h>
#include <float.h>
#include <fstream>
#include <sstream>


#define IPS_URL std::string("https://www.ebi.ac.uk/Tools/services/rest/iprscan5/")


size_t curlWriter(void * buffer, size_t size, size_t nmemb, std::string * s);

std::string submitJob(std::string email, std::string title, std::string sequence);

std::string getJobStatus(std::string jobID);

void getJobResult(std::string jobID, std::string resultType, std::string outFilePath);

std::string getBestMatchTSV(std::string resultTsvFilePath);
