#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "print_info.h"
#include "sra.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_TRIMM_DIR std::string((dl::program_location().parent_path() / fs::path(std::string("../external/Trimmomatic/"))).c_str())
#define PATH_TRIMM std::string((dl::program_location().parent_path() / fs::path(std::string("../external/Trimmomatic/trimmomatic-0.39.jar"))).c_str())
#define TRUSEQ_ADAPT3_PE std::string(PATH_TRIMM_DIR + "/adapters/TruSeq3-PE.fa")
#define TRUSEQ_ADAPT3_SE std::string(PATH_TRIMM_DIR + "/adapters/TruSeq3-SE.fa")
#define TRUSEQ_ALL std::string(PATH_TRIMM_DIR + "/adapters/TruSeq_all.fa")


void run_trimmomatic(std::pair<std::string, std::string> sraRunIn,
                     std::pair<std::string, std::string> sraRunOutP,
                     std::pair<std::string, std::string> sraRunOutU,
                     std::string threads, std::string maxSeedMismatch,
                     std::string minMatchPaired, std::string minMatchAny,
                     std::string windowSize, std::string windowMinQuality,
                     std::string minQualityLead, std::string minQualityTrail,
                     std::string minReadLength, std::string numBpCutFront,
                     bool dispOutput, std::string logFile);
