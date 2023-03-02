#include <iostream>
#include <sstream>
#include <fstream>
#include <string>


class logStream {
  private:
    std::ostream &os1;
    std::ostream &os2;
  public:
    logStream(std::ostream &os1, std::ostream &os2) : os1(os1), os2(os2) {}
    
    template <typename T>
    logStream& operator<<(const T &input) {
      os1 << input;
      os2 << input;
      return *this;
    }
};
