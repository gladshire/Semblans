#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "ini_parse.h"
#include "sra.h"
#include "transcript.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_TRANSD_LONGORFS std::string((dl::program_location().parent_path() / fs::path(std::string("../lib/TransDecoder-TransDecoder-v5.7.0/TransDecoder.LongOrfs"))).c_str())

#define PATH_TRANSD_PREDICT std::string((dl::program_location().parent_path() / fs::path(std::string("../lib/TransDecoder-TransDecoder-v5.7.0/TransDecoder.Predict"))).c_str())

void run_transdecoder(transcript trans, std::string threads);
