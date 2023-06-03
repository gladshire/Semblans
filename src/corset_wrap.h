#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "print_info.h"
#include "ini_parse.h"
#include "sra.h"
#include "transcript.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_CORSET std::string((dl::program_location().parent_path() / fs::path(std::string("../external/corset/corset"))).c_str())

void corset_eq_classes(std::string transPrefix, std::string pathEqClassFile,
                       std::string outDir, bool dispOutput, std::string logFile);
