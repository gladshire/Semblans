#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "ini_parse.h"
#include "print_info.h"
#include "sra.h"
#include "transcript.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_STAR std::string((dl::program_location().parent_path() / fs::path(std::string("../external/star/bin/Linux_x86_64/STAR"))).c_str())

void star_index(std::vector<std::string> fastaFiles, std::string outDir, std::string threads,
                bool dispOutput, std::string logFile);

void star_map(std::string indexPath, std::string outMap,
              std::pair<std::string, std::string> sraReadsIn,
              std::string threads, bool dispOutput, std::string logFile);
