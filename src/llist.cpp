// TODO: getSeq does not properly return sequence -> hash data cannot be updated

#include "llist.h"

linkedList::linkedList() : head(nullptr), tail(nullptr) {}

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

bool linkedList::empty() {
  if (head == nullptr) {
    return true;
  }
  else {
    return false;
  }
}

void linkedList::dump(std::ofstream & outFile) {
  Node * currNode = head;
  std::string currHead;
  std::string currSeq;
  std::string currQual;
  while (currNode != nullptr) {
    currHead = currNode->seqEntry.get_header();
    currSeq = currNode->seqEntry.get_sequence();
    currQual = currNode->seqEntry.get_quality();
    outFile << ">" << currHead << '\n';
    outFile << currSeq << std::endl;
    if (!currQual.empty()) {
      outFile << "+" << '\n';
      outFile << currQual << std::endl;
    }
    currNode = currNode->next;
  }
}

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
