#include "ips_job_man.h"


// IPS sequence annotation job performer
void seqIdJobManager::performJob(std::string email, std::string title, std::string sequence,
                                 std::string annotDir) {
  std::string jobID;
  std::string jobStatus;
  std::string seqBestMatch;
  std::chrono::seconds pollFreq(3);
  // Submit job
  jobID = submitJob(email, title, sequence);
  // Continuously obtain status of job
  while (true) {
    std::this_thread::sleep_for(pollFreq);
    jobStatus = getJobStatus(jobID);
    if (jobStatus == "FINISHED" ||
        jobStatus == "ERROR" ||
        jobStatus == "FAILURE") {
      break;
    }
  }
  if (jobStatus == "FINISHED") {
    // Once finished, obtain GFF3 and TSV result files
    getJobResult(jobID, "gff", annotDir + "/" + title + ".gff3");
    std::this_thread::sleep_for(pollFreq);
    getJobResult(jobID, "tsv", annotDir + "/" + title + ".tsv");
  }
  else {
    // Job failed to complete
  }
}


// Submit sequence to queue of jobs
void seqIdJobManager::submitSeqJob(sequence newSeq) {
  seqJobQueue.push(newSeq);
}


// Initiate job manager on job queue
void seqIdJobManager::performSeqJobs(int numThreads, std::string email, std::string annotDir) {
  std::string title;
  std::string seqData;
  sequence currSeq;

  threadPool annotJobPool;
  if (numThreads > 30) {
    numThreads = 30;
  }

  annotJobPool.start(numThreads);
  while (!seqJobQueue.empty()) {
    currSeq = seqJobQueue.front();
    title = currSeq.get_header();
    seqData = currSeq.get_sequence();
    seqData.pop_back();

    annotJobPool.queueJob([email, title, seqData, annotDir, this]
                          {performJob(email, title, seqData, annotDir);});
    seqJobQueue.pop();
  }
  annotJobPool.stop();

  // Iterate over annotation TSV files in annotation directory, getting best predictions
  fs::directory_iterator fileIter{annotDir};
  std::string currOldSeq;
  std::string currNewSeq;
  while (fileIter != fs::directory_iterator{}) {
    if (fileIter->path().extension() == ".tsv") {
      currOldSeq = std::string(fileIter->path().stem().c_str());
      currNewSeq = getBestMatchTSV(std::string(fileIter->path().c_str()));

      // Emplace (old seq, new seq) pair into map
      newSeqIds.emplace(currOldSeq, currNewSeq);
    }
  }
  std::cout << "Annotation finished. Sequence map should be filled." << std::endl;
}


// Get the map containing the (oldSeqId, newSeqId) pairs
std::map<std::string, std::string> seqIdJobManager::getSeqIds() {
  return newSeqIds;
}



