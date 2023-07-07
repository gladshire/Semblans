#include "print_info.h"

// Given a string, output it to stdout as well as a log file
void logOutput(std::string input, std::string logFile) {
  std::ofstream logStream;
  logStream.open(logFile, std::ios_base::app);
  teedev logger(logStream, std::cout);
  teeStream loggerStream(logger);
  
  loggerStream << input << std::endl;
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
    logOutput("", logFile);
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
  }
}
