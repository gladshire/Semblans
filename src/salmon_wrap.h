#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "ini_parse.h"
#include "sra.h"
#include "transcript.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_SALMON std::string((dl::program_location().parent_path() / fs::path(std::string("../external/salmon-1.9.0/bin/salmon"))).c_str())

void combine_paired(std::vector<SRA> sras, long long int ram_b);

void salmon_index(transcript trans, std::string threads,
                  bool dispOutput, std::string logFile);

void salmon_quant(transcript trans, std::vector<SRA> sras, std::string threads,
                  bool dispOutput, std::string logFile);


