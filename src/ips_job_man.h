#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <queue>
#include <map>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <curl/curl.h>
#include <boost/filesystem.hpp>
#include "ips_client.h"
#include "seq.h"
#include "seq_hash.h"
#include "print_info.h"

class seqJobManager {

  private:
    std::queue<sequence> jobQueue;
    std::map<std::string, sequence> jobsRunning;
    std::map<std::string, sequence> jobsFailed;
    std::map<std::string, std::string> bestMatches;
    int totalJobs;
  public:
    void submitSeqJob(sequence seqJob);
    void startSeqJobs(int numConcurrent, std::string annotDir,
                      std::string email, std::string logFilePath);
    int getNumJobsQueued();
    int getNumJobsRunning();
    int getNumJobsFinished();
    std::map<std::string, std::string> getMatches();
    std::map<std::string, sequence> getJobsFailed();
}; 
