#include "ini_parse.h"
#include "sra.h"


#define PATH_SRA_BIN (std::string)"../lib/sratoolkit.3.0.0-ubuntu64/bin/"
#define PATH_PREFETCH (std::string)PATH_SRA_BIN"prefetch"
#define PATH_FASTERQ (std::string)PATH_SRA_BIN"fasterq-dump"

namespace fs = boost::filesystem;

void prefetch_sra(std::vector<SRA> sras);

void fasterq_sra(std::vector<SRA> sras);

std::vector<SRA> get_sras(const INI_MAP &iniFile);
