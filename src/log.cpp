#include "log.h"

void logOutput(std::string input, std::string logFile) {
  std::ofstream logStream;
  logStream.open(logFile, std::ios_base::app);
  teedev logger(logStream, std::cout);
  teeStream loggerStream(logger);

  loggerStream << input << std::flush;
}
