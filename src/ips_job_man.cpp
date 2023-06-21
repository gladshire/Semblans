#include "ips_job_man.h"


void seqJobManager::submitSeqJob(sequence seqJob) {
  jobQueue.push(seqJob);
}

void seqJobManager::startSeqJobs(int numConcurrent, std::string annotDir,
                                 std::string email, std::string logFilePath) {
  totalJobs = getNumJobsQueued();
  sequence currSeq;
  std::string title;
  std::string seqData;
  std::string jobID;
  std::string jobStatus;
  std::vector<std::string> jobsToRemove;
  std::chrono::seconds pollFreq(5);
  do {
    std::cout << "\nAnnotations queued:   " << getNumJobsQueued() << std::endl;
    std::cout << "Annotations running:  " << getNumJobsRunning() << std::endl;
    std::cout << "Annotations finished: " << getNumJobsFinished() << std::endl;
    
    while (jobsRunning.size() < numConcurrent && !jobQueue.empty()) {
      // Obtain sequence at front of job queue
      currSeq = jobQueue.front();
      title = currSeq.get_header();
      seqData = currSeq.get_sequence();
      seqData.pop_back();

      // Check if outputs for front job already exist
      if (fs::exists(fs::path((annotDir + "/" + title + ".gff3").c_str())) &&
          fs::exists(fs::path((annotDir + "/" + title + ".tsv").c_str()))) {
        std::cout << "Outputs found for: " << title << std::endl;
        continue;
      }

      // Submit sequence as annotation job
      jobID = submitJob(email, title, seqData);

      std::cout << "Submitted: " << jobID << std::endl;

      // Remove sequence from front of job queue
      jobQueue.pop();

      // Add sequence and job ID to list of running jobs
      jobsRunning.emplace(jobID, currSeq);
    }
    for (auto job : jobsRunning) {
      // Get status of running job
      jobStatus = getJobStatus(job.first);

      // If job has ended flag it for removal from list of running jobs
      if (jobStatus != "QUEUED" &&
          jobStatus != "RUNNING") {
        jobsToRemove.push_back(job.first);
        // If job finished successfully, obtain its results
        if (jobStatus == "FINISHED") {
          // Retrieve output in gff3 format
          getJobResult(job.first, "gff", annotDir + "/" +
                       job.second.get_header() + ".gff3");
          // Retrieve output in tsv format
          getJobResult(job.first, "tsv", annotDir + "/" +
                       job.second.get_header() + ".tsv");
        }
        // If job failed, add it to list of failed jobs
        else {
          jobsFailed.emplace(job.first, job.second);
        }
      }
    }
    // Remove flagged jobs from list of running jobs
    for (auto job : jobsToRemove) {
      jobsRunning.erase(job);
    }
    // Wait before submitting further requests
    std::this_thread::sleep_for(pollFreq);
  } while (!jobQueue.empty() && !jobsRunning.empty());
}

int seqJobManager::getNumJobsQueued() {
  return jobQueue.size();
}

int seqJobManager::getNumJobsRunning() {
  return jobsRunning.size();
}

int seqJobManager::getNumJobsFinished() {
  return totalJobs - getNumJobsQueued() - getNumJobsRunning();
}

std::map<std::string, std::string> seqJobManager::getMatches() {
  return bestMatches;
}

std::map<std::string, sequence> seqJobManager::getJobsFailed() {
  return jobsFailed;
}
