#include <boost/filesystem.hpp>
#include "print_info.h"
#include "ini_parse.h"
#include "sra.h"

#define PATH_KRAK2 std::string("../lib/kraken2-2.1.2/kraken2")

std::vector<std::string> get_kraken2_dbs(const INI_MAP &iniFile);

void pre_summary(SRA sra, std::string);

void run_kraken2(std::vector<SRA> sras, std::string threads, std::string db, bool selfPass);

void run_kraken2_dbs(std::vector<SRA> sras, std::string threads, std::vector<std::string> dbs);
