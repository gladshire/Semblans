#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <signal.h>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include <boost/regex.hpp>
#include "log.h"
#include "print_info.h"
#include "ini_parse.h"
#include "sra.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define SEMBLANS_DIR std::string((dl::program_location().parent_path()).c_str()) + "/"


void print_intro(std::string logFile);

void print_help_base();

void print_help_postprocess();
