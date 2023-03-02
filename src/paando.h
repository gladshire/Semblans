#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <sys/ioctl.h>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "ini_parse.h"
#include "sra.h"
#include "logstream.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PAANDO_DIR std::string((dl::program_location().parent_path()).c_str()) + "/"


void print_intro();

void print_help();
