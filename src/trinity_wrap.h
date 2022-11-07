#include <fstream>
#include <boost/filesystem.hpp>
#include "print_info.h"
#include "ini_parse.h"
#include "sra.h"
#include "transcript.h"

#define PATH_TRINITY std::string("../lib/trinityrnaseq-v2.14.0/Trinity")

// Get vector of SRAs to combine into single Trinity assembly
std::vector<SRA> get_sra_to_combine(std::vector<SRA> sras, std::string org_name);

// Concatenate read files for assembly of multiple SRA runs
std::string combine_reads(std::vector<SRA> sras_comb, long long int ram_b);

// Run Trinity for assembly of single SRA run
transcript run_trinity(SRA sra, std::string threads, std::string ram_gb);

// Run Trinity for assembly of several SRA runs
transcript run_trinity_comb(std::vector<SRA> sras_comb, std::string threads, std::string ram_gb);

// Run Trinity for all SRA accessions iteratively
std::vector<transcript> run_trinity_bulk(std::vector<SRA> sras,
                                         std::string threads, std::string ram_gb, bool mult_sra);
