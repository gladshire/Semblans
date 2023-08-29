#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include <boost/regex.hpp>
#include "ini_parse.h"
#include "print_info.h"
#include "sra.h"
#include "transcript.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_SALMON std::string((dl::program_location().parent_path() / fs::path(std::string("../external/salmon/bin/salmon"))).c_str())

void prepareIndexFiles(std::vector<std::string> desiredFiles,
                       std::vector<std::string> decoyFiles,
                       std::string preIndexOutput,
                       std::string decoyFileOutput);

void salmon_index(std::string transIn, std::string transIndex, std::string decoys,
                  std::string threads, bool dispOutput, std::string logFile);

bool runPaired(std::vector<std::pair<std::string, std::string>> sraRunsIn);

void salmon_quant(std::string transIn, std::string transIndex, std::string transQuant,
                  std::vector<std::pair<std::string, std::string>> sraRunsIn,
                  std::string threads, bool dispOutput, std::string logFile);


