#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include <set>
#include "ini_parse.h"
#include "sra.h"
#include "transcript.h"
#include "seq_hash.h"
#include "ncbi_blast.h"
#include "diamond.h"

namespace fs = boost::filesystem;
namespace dl = boost::dll;

#define PATH_TRANSD_LONGORFS std::string((dl::program_location().parent_path() / fs::path(std::string("../external/TransDecoder-5.7.0/TransDecoder.LongOrfs"))).c_str())

#define PATH_TRANSD_PREDICT std::string((dl::program_location().parent_path() / fs::path(std::string("../external/TransDecoder-5.7.0/TransDecoder.Predict"))).c_str())

bool fasta_ok(std::string fastaFile, uintmax_t ram_b);

bool blastpout_ok(std::string blastpFile);

void run_transdecoder(transcript trans, std::string threads, uintmax_t ram_b,
                      std::string dbPath, std::string outDir);
