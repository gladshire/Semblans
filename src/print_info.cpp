#include "print_info.h"


void progressAnim(std::string precedeString, std::string logFile) {

  const std::string anim[] = {".  ", ".. ", "..."};
  int animIndex = 0;

  logOutput(precedeString, logFile);
  while (procRunning) {
    std::cout << "\r" << precedeString;
    std::cout << anim[animIndex] << std::flush;
    animIndex = (animIndex + 1) % 3;
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  }
  std::cout << "\r" << precedeString;
  std::cout << "    " << std::endl;
}


// Given an SRA run object, output basic information to stdout and
// a log file
void summarize_sing_sra(SRA sra, std::string logFile, int margin) {
  std::ofstream logStream;
  logStream.open(logFile, std::ios_base::app);
  teedev logger(logStream, std::cout);
  teeStream loggerStream(logger);
  std::string filePrefix1 = sra.get_file_prefix().first;
  std::string filePrefix2 = sra.get_file_prefix().second;
  std::string margStr;
  for (int i = 0; i < margin; i++) {
    margStr += " ";
  }
  if (sra.is_paired()) {
    loggerStream << margStr << "Paired-end run:" << std::endl;
    loggerStream << margStr << filePrefix1 << std::endl;
    loggerStream << margStr << filePrefix2 << std::endl;
  }
  else {
    loggerStream << margStr << "Single-end run:" << std::endl;
    loggerStream << margStr << filePrefix1 << std::endl;
  }
  logStream.close();
}

// Given a vector of SRA run objects, briefly summarize basic information for
// each in stdout and a log file
void summarize_all_sras(const std::vector<SRA> & sras, std::string logFile,
                        int margin) {
  for (auto sra : sras) {
    summarize_sing_sra(sra, logFile, margin);
    logOutput("\n", logFile);
  }
}

// Given a transcript object, output basic information to stdout and
// a log file
void summarize_sing_trans(transcript trans, std::string logFile, int margin) {
  std::ofstream logStream;
  logStream.open(logFile, std::ios_base::app);
  teedev logger(logStream, std::cout);
  teeStream loggerStream(logger);
  std::string filePrefix = trans.get_file_prefix();
  std::string margStr;
  for (int i = 0; i < margin; i++) {
    margStr += " ";
  }
  loggerStream << margStr << "Transcript:" << std::endl;
  loggerStream << margStr << filePrefix << std::endl;
  logStream.close();
}

// Given a vector of transcript objects, briefly summarize basic information for
// each in stdout and a log file
void summarize_all_trans(const std::vector<transcript> & transVec, std::string logFile,
                         int margin) {
  for (auto trans : transVec) {
    summarize_sing_trans(trans, logFile, margin);
    logOutput("\n", logFile);
  }
}
