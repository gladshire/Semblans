#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "print_info.h"
#include "sra.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_TRIMM_DIR std::string((dl::program_location().parent_path() / fs::path(std::string("../lib/Trimmomatic-0.39/"))).c_str())
#define PATH_TRIMM std::string((dl::program_location().parent_path() / fs::path(std::string("../lib/Trimmomatic-0.39/trimmomatic-0.39.jar"))).c_str())
#define TRUSEQ_ADAPT3_PE std::string(PATH_TRIMM_DIR + "/adapters/TruSeq3-PE.fa")
#define TRUSEQ_ADAPT3_SE std::string(PATH_TRIMM_DIR + "/adapters/TruSeq3-SE.fa")
#define TRUSEQ_ALL std::string(PATH_TRIMM_DIR + "/adapters/TruSeq_all.fa")


void run_trimmomatic(std::vector<SRA> sras, std::string threads);
