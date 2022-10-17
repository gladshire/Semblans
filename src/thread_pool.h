#include <thread>
#include <mutex>
#include <functional>
#include <vector>
#include <queue>
#include <condition_variable>

class threadPool {
  public:
    void start(int threadNum);
    void queueJob(const std::function<void()>& job);
    void stop();
    bool busy();

  private:
    void threadLoop();

    bool endProc = false;
    std::mutex queue_mutex;
    std::condition_variable mutex_condition;
    std::vector<std::thread> threads;
    std::queue<std::function<void()>> jobs;
};
