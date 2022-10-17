#include <boost/filesystem.hpp>
#include "print_info.h"
#include "sra.h"

#define PATH_TRIMM std::string("../lib/Trimmomatic-0.39/trimmomatic-0.39.jar")
#define TRUSEQ_ADAPT3_PE std::string("../lib/Trimmomatic-0.39/adapters/TruSeq3-PE.fa")
#define TRUSEQ_ADAPT3_SE std::string("../lib/Trimmomatic-0.39/adapters/TruSeq3-SE.fa")
#define TRUSEQ_ALL std::string("../lib/Trimmomatic-0.39/adapters/TruSeq_all.fa")

namespace fs = boost::filesystem;

void run_trimmomatic(std::vector<SRA> sras, std::string threads);
