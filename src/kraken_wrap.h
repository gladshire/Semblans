#include <boost/filesystem.hpp>
#include "sra.h"

#define PATH_KRAK2 std::string("../lib/kraken2-2.1.2/kraken2")

void run_kraken2(std::vector<SRA> sras, std::string threads, std::string db, bool selfPass);

void run_kraken2_dbs(std::vector<SRA> sras, std::string threads, std::vector<std::string> dbs);
