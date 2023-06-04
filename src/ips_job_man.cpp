#include "ips_job_man.h"


// IPS sequence annotation job performer
void seqIdJobManager::performJob(std::string email, std::string title, std::string sequence) {
  std::string jobID;
  std::string jobStatus;
  std::string seqBestMatch;
  // Submit job
  jobID = submitJob(email, title, sequence);

  // Continuously obtain status of job
  while (true) {
    jobStatus = getJobStatus(jobID);
    if (jobStatus == "FINISHED" ||
        jobStatus == "ERROR" ||
        jobStatus == "FAILURE") {
      break;
    }
  }
  if (jobStatus == "FINISHED") {
    // Once finished, obtain TSV result file
    getJobResult(jobID, "tsv", title + ".tsv");
    // Get best sequence match
    seqBestMatch = getBestMatchTSV(title + ".txv");
    // Emplace new (oldSeqId, newSeqId) entry into map
    newSeqIds.emplace(title, seqBestMatch);
  }
}


// Submit sequence to queue of jobs
void seqIdJobManager::submitSeqJob(sequence newSeq) {
  seqJobQueue.push(newSeq);
}


// Initiate job manager on job queue
void seqIdJobManager::startSeqJobs(int numThreads, std::string email) {
  std::string title;
  std::string seqData;
  sequence currSeq;
  while (!seqJobQueue.empty()) {
    currSeq = seqJobQueue.front();
    title = currSeq.get_header();
    seqData = currSeq.get_sequence();
    seqJobPool->queueJob([&] {performJob(email, title, seqData);});
    seqJobQueue.pop();
  }
  if (numThreads > 30) {
    numThreads = 30;
  }
  seqJobPool->start(numThreads);
  while (seqJobPool->busy()) {
    std::cout << "Thread pool busy" << std::endl;
  }
  std::cout << "Done. Seq map should be filled" << std::endl;
}


// Get the map containing the (oldSeqId, newSeqId) pairs
std::map<std::string, std::string> seqIdJobManager::getSeqIds() {
  return newSeqIds;
}



