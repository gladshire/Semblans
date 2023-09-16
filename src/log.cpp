#include "log.h"

// Utility function to replace one character with another in a text file
void replaceChar(std::string inFilePath, char oldChar, char newChar) {
  std::ifstream inLogFile(inFilePath);
  std::ofstream outLogFile(inFilePath + ".tmp");
   
  char currChar;
  while (inLogFile.get(currChar)) {
    if (currChar != oldChar) {
      outLogFile.put(currChar);
    }
    else {
      outLogFile.put(newChar);
    }
  }

  inLogFile.close();
  outLogFile.close();

  std::remove(inFilePath.c_str());
  std::rename((inFilePath + ".tmp").c_str(), inFilePath.c_str());
}

// Primary output function
// Given a string, it prints it to std::out and writes it to a log file
void logOutput(std::string input, std::string logFile) {
  std::ofstream logStream;
  logStream.open(logFile, std::ios_base::app);
  teedev logger(logStream, std::cout);
  teeStream loggerStream(logger);

  loggerStream << input << std::flush;
}
