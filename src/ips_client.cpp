#include "ips_client.h"


size_t curlWriter(void * buffer, size_t size, size_t nmemb, std::string * s) {
  size_t newLen = size * nmemb;
  try {
    s->append((char*)buffer, newLen);
  }
  catch (std::bad_alloc &e) {
    // memory problem
    return 0;
  }
  return newLen;
}

// Function to submit annotation job to interproscan rest web API
std::string submitJob(std::string email, std::string title, std::string sequence) {
  std::string jobID;
  std::string postData;
  std::string runURL = IPS_URL + "run/";

  CURL * curl;
  CURLcode result;

  // Initialize libcurl environment
  curl_global_init(CURL_GLOBAL_ALL);

  // Instantiate curl handle
  curl = curl_easy_init();
  if (curl) {
    // Set URL to interproscan API
    curl_easy_setopt(curl, CURLOPT_URL, runURL.c_str());

    // Set request to POST
    curl_easy_setopt(curl, CURLOPT_POST, 1L);

    // Construct post data
    postData = "email=" + email + "&" +
               "title=" + title + "&" +
               "sequence=" + sequence;
    
    // Set size of post data
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, postData.length());

    // Set request data
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, postData.c_str());

    // Set callback function to handle response
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curlWriter);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &jobID);

    // Execute post request
    result = curl_easy_perform(curl);

    // Check for errors
    if (result != CURLE_OK) {
      //logOutput("curl_easy_perform() failed: " + curl_easy_strerror(result), logFile);
      std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(result) << std::endl;
    }
    curl_easy_cleanup(curl);
  }
  curl_global_cleanup();
  return jobID;
}


// Retrieve current status of job
//
//   Possible status values:
//
//     QUEUED:    Job is queued
//     RUNNING:   Job being processed
//     FINISHED:  Job is complete
//     ERROR:     Job experienced error 
//     FAILURE:   Job failed to complete
//     NOT_FOUND: Job cannot be found (invalid job id)

std::string getJobStatus(std::string jobID) {
  std::string jobStatus;
  std::string statusURL = IPS_URL + "status/" + jobID;

  CURL * curl;
  CURLcode result;

  // Initialize libcurl environment
  curl_global_init(CURL_GLOBAL_ALL);

  // Instantiate curl handle
  curl = curl_easy_init();

  if (curl) {
    // Set URL to interproscan API
    curl_easy_setopt(curl, CURLOPT_URL, statusURL.c_str());

    // Set callback function to handle response
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curlWriter);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &jobStatus);

    // Execute request
    result = curl_easy_perform(curl);

    // Check for errors
    if (result != CURLE_OK) {
      std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(result) << std::endl;
    }
    curl_easy_cleanup(curl);
  }
  curl_global_cleanup();
  return jobStatus;
}


// Retrieve result of job
void getJobResult(std::string jobID, std::string resultType, std::string outFilePath) {
  std::ofstream outFile;
  std::string jobResult;
  std::string resultURL = IPS_URL + "result/" + jobID + "/" + resultType;

  CURL * curl;
  CURLcode result;

  // Initialize libcurl environment
  curl_global_init(CURL_GLOBAL_ALL);

  // Instantiate curl handle
  curl = curl_easy_init();

  if (curl) {
    // Set URL to interproscan API
    curl_easy_setopt(curl, CURLOPT_URL, resultURL.c_str());

    // Set callback function to handle response
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curlWriter);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &jobResult);

    // Execute request
    result = curl_easy_perform(curl);

    // Check for errors
    if (result != CURLE_OK) {
      std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(result) << std::endl;
    }
    curl_easy_cleanup(curl);
  }
  curl_global_cleanup();
  outFile.open(outFilePath);
  outFile << jobResult;
}


// Obtain likeliest ID from interproscan result TSV file
std::string getBestMatchTSV(std::string resultTsvFilePath) {
  std::ifstream resultTsvFile(resultTsvFilePath);
  std::string currLine;
  std::string currTok;
  std::string eValStr;
  std::string matchId;
  std::string matchIdBest;
  double eValMin = DBL_MAX;
  double eValCurr = DBL_MAX;
  size_t tabPos;
  int linNum = 0;
  while (getline(resultTsvFile, currLine)) {
    // Split line by tab
    //std::cout << currLine << std::endl;
    int colNum = 0;
    while ((tabPos = currLine.find("\t")) != std::string::npos) {
      currTok = currLine.substr(0, tabPos);
      if (colNum == 5) {
        matchId = currTok;
      }
      if (colNum == 8) {
        eValStr = currTok;
      }
      currLine.erase(0, tabPos + 1);
      colNum++;
    }
    std::istringstream os(eValStr);
    if (eValStr != "-" && eValStr != "0.0") {
      os >> eValCurr;
    }
    // Want first number number found to be set as min
    // Any subsequent number less than last set as min
    if (linNum == 0 || (eValCurr < eValMin && eValStr != "-")) {
      eValMin = eValCurr;
      matchIdBest = matchId;
    }
    linNum++;
  }
  return matchIdBest;
}
