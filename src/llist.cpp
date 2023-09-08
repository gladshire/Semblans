// Implementation file for simple, singly-linked list class containing sequence
// data

#include "llist.h"

// Default constructor for linked list
linkedList::linkedList() : head(nullptr), tail(nullptr) {}

// Insert function for appending a sequence node to the linked list
void linkedList::insert(const sequence & seqData) {
  Node * newNode = new Node(seqData);
  Node * currNode = head;
  if (head == nullptr) {
    head = newNode;
  }
  else {
    while (currNode->next != nullptr) {
      currNode = currNode->next;
    }
    currNode->next = newNode;
  }
}

// Remove function for deleting a sequence node with a given header from
// the linked list
bool linkedList::remove(const std::string & header) {
  Node * currNode = head;
  Node * lastNode = nullptr;
  std::string currHead;
  while (currNode != nullptr) {
    currHead = currNode->seqEntry.get_header();
    if (currHead.substr(0, currHead.find(" ")) == header ||
        currHead == header) {
      if (lastNode != nullptr) {
        lastNode->next = currNode->next;
      }
      else {
        head = currNode->next;
      }
      delete currNode;
      return true;
    }
    lastNode = currNode;
    currNode = currNode->next;
  }
  return false;
}

// Getter function for returning a sequence with a given header from the
// linked list
sequence linkedList::getSeq(const std::string & header) {
  Node * currNode = head;
  std::string currHead;
  while (currNode != nullptr) {
    currHead = currNode->seqEntry.get_header();
    if (currHead.substr(0, currHead.find(" ")) == header ||
        currHead == header) {
      return sequence(currNode->seqEntry.get_header(),
                      currNode->seqEntry.get_sequence(),
                      currNode->seqEntry.get_quality());
    }
    currNode = currNode->next;
  }
  return sequence();
}

// Setter function for updating the header of a sequence node with
// a given 'old header' to a 'new header' in the linked list
void linkedList::setSeqHead(const std::string & oldHead,
                            const std::string & newHead) {
  Node * currNode = head;
  while (currNode != nullptr) {
    if (currNode->seqEntry.get_header() == oldHead) {
      currNode->seqEntry.set_header(newHead);
      return;
    }
    currNode = currNode->next;
  }
}

// Utility function for updating all headers in the linked list based on
// a given prefix
// Used primarily for updating FASTA transcript headers post-assembly
void linkedList::updateSeqHead(const std::string & newPrefix) {
  Node * currNode = head;
  std::string currHead;
  std::string newHead;
  while (currNode != nullptr) {
    currHead = currNode->seqEntry.get_header();
    newHead = newPrefix + currHead.substr(currHead.find("_"));
    setSeqHead(currHead, newHead);
    currNode = currNode->next;
  }
}

// Logical function for determining whether a sequence node with a given
// header exists in the linked list
bool linkedList::exists(const std::string & header) {
  Node * currNode = head;
  std::string currHead;
  while (currNode != nullptr) {
    currHead = currNode->seqEntry.get_header();
    if (currHead.substr(0, currHead.find(" ")) == header ||
        currHead == header) {
      return true;
    }
    currNode = currNode->next;
  }
  return false;
}

// Logical function for determining whether the linked list is empty
bool linkedList::empty() {
  if (head == nullptr) {
    return true;
  }
  else {
    return false;
  }
}

// Function for dumping all sequence data in the linked list to a
// specified output stream
void linkedList::dump(std::ofstream & outFile) {
  Node * currNode = head;
  std::string currHead;
  std::string currSeq;
  std::string currQual;
  while (currNode != nullptr) {
    currHead = currNode->seqEntry.get_header();
    currSeq = currNode->seqEntry.get_sequence();
    currQual = currNode->seqEntry.get_quality();
    if (!currQual.empty()) {
      outFile << "@";
    }
    else {
      outFile << '>';
    }
    outFile << currHead << '\n';
    outFile << currSeq << std::endl;
    if (!currQual.empty()) {
      outFile << "+" << '\n';
      outFile << currQual << std::endl;
    }
    currNode = currNode->next;
  }
}

// Function for deleting and freeing all sequence data contained in the
// linked list
void linkedList::clear() {
  Node * tmpNode = head;
  Node * currNode = head;
  while (currNode != nullptr) {
    tmpNode = currNode;
    currNode = currNode->next;
    delete tmpNode;
  }
  head = nullptr;
  tail = nullptr;
}
