#include "log.h"



// Utility function to replace one character with another in a text file
void replaceChar(std::string inFilePath, char oldChar, char newChar) {
  std::ifstream inLogFile(inFilePath);
  std::ofstream outLogFile(inFilePath + ".tmp");
   
  char currChar;
  char lastChar;
  while (inLogFile.get(currChar)) {
    if (currChar == oldChar) {
      if (currChar != lastChar) {
        outLogFile.put(newChar);
      }
    }
    else {
      outLogFile.put(currChar);
    }
    lastChar = currChar;
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
  char timeBuffer[100];
  time_t currTime;
  logStream.open(logFile, std::ios_base::app);
  teedev logger(logStream, std::cout);
  teeStream loggerStream(logger);

  time(&currTime);
  strftime(timeBuffer, 100, "[%R:%S %F]  ", localtime(&currTime));

  //loggerStream << timeBuffer << std::flush;

  loggerStream << input << std::flush;
}

void printVertEllipse(std::string logFile, int numLines) {
  for (int i = 0; i < numLines; i++) {
    logOutput("                            .                           \n", logFile);
  }
}

void printBreakLine(std::string logFile, int leadingSpace, int length) {
  for (int i = 0; i < leadingSpace; i++) {
    logOutput(" ", logFile);
  }
  for (int i = 0; i < length; i++) {
    logOutput("─", logFile);
  }
  logOutput("\n", logFile);
}

// Return string of percentage to a given precision
std::string getPercent(float valPercent, int precision) {
  std::stringstream percentStream;
  percentStream << std::fixed << std::setprecision(precision) << valPercent;
  return percentStream.str();
}

// Print error code of a program call and then exit
void checkExitSignal(int commandResult, std::string logFile) {
  if (WIFSIGNALED(commandResult)) {
    system("setterm -cursor on");
    logOutput("\nExited with signal " + std::to_string(WTERMSIG(commandResult)) + "\n", logFile);
    exit(1);
  }
}

// Given a file path, return its prefix with no extensions
std::string removeExtensions(std::string filePath) {
  fs::path filePrefix(fs::path(filePath.c_str()));
  while (!filePrefix.extension().empty()) {
    filePrefix = filePrefix.stem();
  }
  std::string fileStem(filePrefix.c_str());
  return fileStem;
}
