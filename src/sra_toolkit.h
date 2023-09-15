#pragma once
#include "sra.h"
#include "print_info.h"
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>

namespace dl = boost::dll;
namespace fs = boost::filesystem;

#define PATH_SRA_BIN std::string((dl::program_location().parent_path() / fs::path(std::string("../external/sratoolkit/bin/"))).c_str())
#define PATH_PREFETCH std::string(PATH_SRA_BIN + "prefetch")
#define PATH_FASTERQ std::string(PATH_SRA_BIN + "fasterq-dump")
#define PATH_PIGZ std::string((dl::program_location().parent_path() / fs::path(std::string("../external/pigz/pigz"))).c_str())

std::vector<SRA> get_sras(const INI_MAP &iniFile, bool dispOutput, bool compressFiles,
                          std::string logFile);

void prefetch_sra(SRA sra, bool dispOutput,
                  std::string logFile);

void fasterq_sra(SRA sra, std::string threads,
                 bool dispOutput, bool compressOutput,
                 std::string logFile);

void align_file_buffer(std::istream & inStream1, std::istream & inStream2,
                       char * inStream1Data, char * inStream2Data,
                       std::streamsize & s1, std::streamsize & s2);
