#include "seq.h"

// Default constructor for sequence object
sequence::sequence() {
  header = "";
  sequenceData = "";
  quality = "";
  int numBp = -1;
}

// Constructor for sequence object allowing definition of its header and
// sequence data
sequence::sequence(std::string header, std::string sequenceData) {
  this->header = header;
  size_t nlPos = sequenceData.find("\n");
  this->sequenceData = sequenceData;
  this->quality = "";
  this->id = "";
  numBp = sequenceData.length();
}

// Constructor for sequence object allowing definition of its header, sequence data, and
// quality data
sequence::sequence(std::string header, std::string sequenceData, std::string quality) {
  this->header = header;
  size_t nlPos = sequenceData.find('\n');
  this->sequenceData = sequenceData;
  this->quality = quality;
  this->id = "";
  numBp = sequenceData.length();
}

// Copy constructor for sequence object
sequence::sequence(const sequence & seq) {
  header = seq.header;
  sequenceData = seq.sequenceData;
  quality = seq.quality;
}

// Return sequence object's header
std::string sequence::get_header() {
  //std::cout << header << std::endl;
  return this->header;
}

// Return sequence object's sequence data
std::string sequence::get_sequence() {
  return sequenceData;
}

// Return sequence object's quality data
std::string sequence::get_quality() {
  return quality;
}

// Return sequence object's ID
std::string sequence::get_id() {
  return id;
}

// Set new header for sequence object
void sequence::set_header(std::string newHeader) {
  header = newHeader;
}

// Set new ID for sequence object
void sequence::set_id(std::string newId) {
  id = newId;
}
