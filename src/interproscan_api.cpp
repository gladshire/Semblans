//#include "interproscan_api.h"
#include <string>
#include <vector>
#include <iostream>
#include <cstring>
#include <curl/curl.h>


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


std::string submitJob(std::string email, std::string title, std::string sequence) {
  std::string jobID;
  std::string postData;
  std::string responseStr;

  CURL * curl;
  CURLcode result;

  // Initialize libcurl environment
  curl_global_init(CURL_GLOBAL_ALL);

  // Instantiate curl handle
  curl = curl_easy_init();

  if (curl) {
    // Set URL to interproscan API
    curl_easy_setopt(curl, CURLOPT_URL, IPS_URL);

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
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &responseStr);

    // Execute post request
    result = curl_easy_perform(curl);

    // Check for errors
    if (result != CURLE_OK) {
      //logOutput("curl_easy_perform() failed: " + curl_easy_strerror(result), logFile);
      std::cout << "curl_easy_perform() failed: " + curl_easy_strerror(result) << std::endl;
    }
    curl_easy_cleanup(curl);
  }
  curl_global_cleanup();

  return jobID;
}


int main() {
  std::string email = "gladshire@gmail.com";
  std::string title = "test";
  std::string sequence = ">tr|B4Q421|B4Q421_DROSI Phospholipase B1, membrane-associated OS=Drosophila simulans OX=7240 GN=Dsim\\GD23360 PE=3 SV=1\n
                          MASITYALCLCSLFILFSPLVSSNRRVRRQNRSDRLLADIGVRAQSYDPRVPENGIQQYT\n
                          DIDQDLRHLFLNTRQTTLRWALNNIDALSSRGRREGKLQVPVSKKVPFLCPTNDTRSPSP\n
                          PTSIEHLRPGDIDIIAAFGDSLSAGNGILSNNAIDMINEFRGLTFSGGGLGNWRRFVTLP\n
                          NILKIFNPKLYGFAVSNSLVINHRSSRLNIAEPMIMSRDLPFQARVLIDLLRRDRHVDMK\n
                          RHWKLLTVYVGNNDICSDLCHWDTPQSFLDQHARDLRQAFRLLRDHVPRLLINLIVVPNI\n
                          PLVLSTMTQVPLQCFVVHRVGCHCLINDRLNRTQFNERMDTLIRWQQLDMEIARLPEFHR\n
                          QDFAIVAHPMLTKVTAPTLPDGSTDWRFFSHDCFHFSQRGHAIISNLLWNSMLLPDDQKP\n
                          RPSVVPELFERVVCPTAEQPYFVVRPS\n";

  std::string jobID;
  jobID = submitJob(email, title, sequence);

}
