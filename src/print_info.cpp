#include "print_info.h"

void logOutput(std::string input, std::string logFile) {
  std::ofstream logStream;
  logStream.open(logFile, std::ios_base::app);
  teedev logger(logStream, std::cout);
  teeStream loggerStream(logger);
  
  loggerStream << input << std::endl;
}

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
    loggerStream << margStr << filePrefix2 << "\n" << std::endl;
  }
  else {
    loggerStream << margStr << "Single-end run:" << std::endl;
    loggerStream << margStr << filePrefix1 << "\n" << std::endl;
  }
  logStream.close();
}

void summarize_all_sras(const std::vector<SRA> & sras, std::string logFile,
                        int margin) {
  for (auto sra : sras) {
    summarize_sing_sra(sra, logFile, margin);
  }
}

