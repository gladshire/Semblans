#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <queue>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <curl/curl.h>
#include <boost/filesystem.hpp>
#include "ips_client.h"
#include "seq.h"
#include "seq_hash.h"
#include "thread_pool.h"


// Job manager structure for interproscan annotation jobs
//
//   Functional goals:
//     - Receive sequence objects for annotation
//     - Organize objects into queue of jobs
//     - Run a certain number of jobs in parallel (max 30)
//     - Maintain track of job statuses
//

//   Implementation:
//     - Place sequence jobs into queue
//     - Queue jobs with thread pool via job-performing function
//       - When a job complete, emplace (oldSeqHeader, newSeqHeader) into newSeqIds map
//     - When jobs finished, obtain newSeqIds map

class seqIdJobManager {
  private:
    std::queue<sequence> seqJobQueue;
    threadPool * seqJobPool;
    std::map<std::string, std::string> newSeqIds;
    void performJob(std::string email, std::string title, std::string sequence);
  public:
    void submitSeqJob(sequence newSeq);
    void startSeqJobs(int numThreads, std::string email);
    std::map<std::string, std::string> getSeqIds();
};
