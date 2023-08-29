// Header file for the simple, singly-linked list class for storing
// FASTA or FASTQ sequence data

#pragma once
#include <iostream>
#include <fstream>
#include "seq.h"


// Struct definition for a sequence node, capable of storing a FASTA/FASTQ
// entry's HEADER, SEQUENCE, and QUALITY strings
struct Node {
  sequence seqEntry;
  Node * next;

  Node(const sequence & seqData) :
       seqEntry(seqData), next(NULL) {}

  Node(const std::string & header,
       const std::string & seqData,
       const std::string & quality) :
       seqEntry(header, seqData, quality), next(NULL) {}

  Node(const std::string & header,
       const std::string & seqData) :
       seqEntry(header, seqData), next(NULL) {}

};

// Class definition for the singly-linked list for storing FASTA/FASTQ data
class linkedList {
  private:
    Node * head;
    Node * tail;

  public:
    linkedList();
    void insert(const sequence & seqData);
    bool remove(const std::string & header);
    sequence getSeq(const std::string & header);
    void setSeqHead(const std::string & oldHead,
                    const std::string & newHead);
    void updateSeqHead(const std::string & newPrefix);
    bool exists(const std::string & header);
    bool empty();
    void dump(std::ofstream & outFile);
    void clear();
};


