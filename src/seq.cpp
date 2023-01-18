#include "seq.h"

sequence::sequence() {
  header = "";
  sequenceData = "";
  quality = "";
  int numBp = -1;
  int numLine = -1;
}

sequence::sequence(std::string header, std::string sequenceData) {
  this->header = header;
  size_t nlPos = sequenceData.find('\n');
  bpPerLine = nlPos + 1;
  numLine = 0;
  while (nlPos != std::string::npos) {
    numLine++;
    sequenceData.erase(nlPos, 1);
    nlPos = sequenceData.find('\n');
  }
  this->sequenceData = sequenceData;
  numBp = sequenceData.length();
}

std::string sequence::get_header() {
  return this->header;
}

std::string sequence::get_sequence() {
  return this->sequenceData;
}

