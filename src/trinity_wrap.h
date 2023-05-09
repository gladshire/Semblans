#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "print_info.h"
#include "ini_parse.h"
#include "sra.h"
#include "transcript.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_TRINITY std::string((dl::program_location().parent_path() / fs::path(std::string("../external/trinityrnaseq-2.15.1/Trinity"))).c_str())

// Get vector of SRAs to combine into single Trinity assembly
std::vector<SRA> get_sra_to_combine(std::vector<SRA> sras, std::string org_name);

// Concatenate read files for assembly of multiple SRA runs
std::string combine_reads(std::vector<std::pair<std::string, std::string>> sraRuns,
                          std::string outFile, long long int ram_b, std::string logFile);

// Run Trinity for assembly of single SRA run
void run_trinity(std::pair<std::string, std::string> sraRun, std::string outFile,
                 std::string threads, std::string ram_gb,
                 bool dispOutput, std::string logFile);

// Run Trinity for assembly of several SRA runs
void run_trinity_comb(std::vector<std::pair<std::string, std::string>> sraRuns,
                      std::string outFile,
                      std::string threads, std::string ram_gb,
                      bool dispOutput, std::string logFile);
