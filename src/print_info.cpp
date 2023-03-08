#include "print_info.h"

void logOutput(std::string input, std::string logFile) {
  std::ofstream logStream;
  logStream.open(logFile, std::ios_base::app);
  teedev logger(logStream, std::cout);
  teeStream loggerStream(logger);
  
  loggerStream << input << std::endl;
}

void summarize_sing_sra(SRA sra, std::string logFile) {
  std::ofstream logStream;
  logStream.open(logFile, std::ios_base::app);
  teedev logger(logStream, std::cout);
  teeStream loggerStream(logger);
  std::string filePrefix1 = sra.get_file_prefix().first;
  std::string filePrefix2 = sra.get_file_prefix().second;
  if (sra.is_paired()) {
    loggerStream << "  Paired-end run:" << std::endl;
    loggerStream << "  " << filePrefix1 << std::endl;
    loggerStream << "  " << filePrefix2 << "\n" << std::endl;
  }
  else {
    loggerStream << "  Single-end run:" << std::endl;
    loggerStream << "  " << filePrefix1 << "\n" << std::endl;
  }
  logStream.close();
}

void summarize_all_sras(const std::vector<SRA> & sras, std::string logFile) {
  for (auto sra : sras) {
    summarize_sing_sra(sra, logFile);
  }
}

