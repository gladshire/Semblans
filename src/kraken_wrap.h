#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "print_info.h"
#include "ini_parse.h"
#include "sra.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_KRAK2 std::string((dl::program_location().parent_path() / fs::path(std::string("../external/kraken2-2.1.2/kraken2"))).c_str())

std::vector<std::string> get_kraken2_dbs(const INI_MAP &iniFile);

std::string get_kraken2_conf(const INI_MAP &iniFile);

void pre_summary(SRA sra, std::string);

void run_kraken2(const std::vector<SRA> & sras, std::string threads, std::string db,
                 std::string conf_threshold, bool selfPass, bool dispOutput, std::string logFile);

void run_kraken2_dbs(const std::vector<SRA> & sras, std::string threads, std::vector<std::string> dbs,
                     std::string conf_threshold, bool dispOutput, std::string logFile);
