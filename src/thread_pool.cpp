#include "thread_pool.h"

void threadPool::start(int threadNum) {
  threads.resize(threadNum);
  for (int i = 0; i < threadNum; i++) {
    threads.at(i) = std::thread(&threadPool::threadLoop, this);
  }
}


void threadPool::threadLoop() {
  while (true) {
    std::function<void()> job;
    {
      std::unique_lock<std::mutex> lock(queue_mutex);
      mutex_condition.wait(lock, [this] {
                           return (!jobs.empty() || endProc);
                           });
      if (jobs.empty() && endProc) {
        return;
      }
      job = jobs.front();
      jobs.pop();
    }
    job();
  }
}


void threadPool::queueJob(const std::function<void()>& job) {
  {
    std::unique_lock<std::mutex> lock(queue_mutex);
    jobs.push(job);
  }
  mutex_condition.notify_one();
}


bool threadPool::busy() {
  bool poolBusy;
  {
    std::unique_lock<std::mutex> lock(queue_mutex);
    poolBusy = jobs.empty();
  }
  return poolBusy;
}

void threadPool::stop() {
  {
    std::unique_lock<std::mutex> lock(queue_mutex);
    endProc = true;
  }
  mutex_condition.notify_all();
  for (std::thread& t : threads) {
    t.join();
  }
  threads.clear();
}
