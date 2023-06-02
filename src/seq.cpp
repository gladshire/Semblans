#include "seq.h"

sequence::sequence() {
  header = "";
  sequenceData = "";
  quality = "";
  int numBp = -1;
}

sequence::sequence(std::string header, std::string sequenceData) {
  this->header = header;
  size_t nlPos = sequenceData.find('\n');
  this->sequenceData = sequenceData;
  this->id = "";
  numBp = sequenceData.length();
}

sequence::sequence(const sequence & seq) {
  header = seq.header;
  sequenceData = seq.sequenceData;
  quality = seq.quality;
  numBp = seq.numBp;
  id = seq.id;
}

std::string sequence::get_header() {
  return this->header;
}

std::string sequence::get_sequence() {
  return this->sequenceData;
}

std::string sequence::get_id() {
  return this->id;
}
