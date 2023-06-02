#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "ini_parse.h"
#include "print_info.h"
#include "sra.h"
#include "transcript.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_SALMON std::string((dl::program_location().parent_path() / fs::path(std::string("../external/salmon-1.10.1/bin/salmon"))).c_str())

void combine_paired(std::vector<SRA> sras, long long int ram_b);

void salmon_index(std::string transIn, std::string transIndex,
                  std::string threads, bool dispOutput, std::string logFile);

bool runPaired(std::vector<std::pair<std::string, std::string>> sraRunsIn);

void salmon_quant(std::string transIn, std::string transIndex, std::string transQuant,
                  std::vector<std::pair<std::string, std::string>> sraRunsIn,
                  std::string threads, bool dispOutput, std::string logFile);


