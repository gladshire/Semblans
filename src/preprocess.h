#pragma once
#include <iostream>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <thread>
#include <chrono>
#include <atomic>
#include <sys/ioctl.h>
#include <unistd.h>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include <boost/regex.hpp>
#include "sra.h"
#include "ini_parse.h"
#include "sra_toolkit.h"
#include "fastqc_wrap.h"
#include "rcorr_wrap.h"
#include "rem_unfixable.h"
#include "trimm_wrap.h"
#include "kraken_wrap.h"
#include "rem_overrep.h"
#include "log.h"
#include "print_info.h"

#define PATH_PIGZ std::string((dl::program_location().parent_path() / fs::path(std::string("../external/pigz/pigz"))).c_str())


std::vector<std::pair<std::string, std::string>> retrieve_sra_data(std::vector<SRA> & sras, std::string threads,
                       bool dispOutput, bool retainInterFiles,
                       bool compressOutput, std::string logFile);

bool stringToBool(std::string boolStr);

void print_help();

void fastqcBulk1(std::vector<SRA> & sras, std::string threads, bool dispOutput,
                 std::string logFilePath);

void fastqcBulk2(std::vector<SRA> & sras, std::string threads, bool dispOutput,
                 std::string logFilePath, std::string outDir, const INI_MAP & cfgIni);

std::vector<std::pair<std::string, std::string>> errorCorrBulk(std::vector<SRA> & sras, std::string threads,
                   bool dispOutput, bool retainInterFiles, bool compressFiles,
                   std::string logFilePath, std::string outDir,
                   const INI_MAP & cfgIni);

std::vector<std::pair<std::string, std::string>> remUnfixBulk(std::vector<SRA> & sras, std::string threads, std::string ram_gb,
                  bool dispOutput, bool retainInterFiles, bool compressFiles,
                  std::string logFilePath, std::string outDir,
                  const INI_MAP & cfgIni);

std::vector<std::pair<std::string, std::string>> trimBulk(std::vector<SRA> & sras, std::string threads,
              bool dispOutput, bool retainInterFiles,
              std::string logFilePath, std::string outDir,
              const INI_MAP & cfgIni);

std::vector<std::pair<std::string, std::string>> filtForeignBulk(std::vector<SRA> & sras, std::vector<std::string> krakenDbs,
                     std::string threads, bool dispOutput, bool compressFiles, bool retainInterFiles,
                     std::string logFilePath, std::string outDir, const INI_MAP & cfgIni);

std::vector<std::pair<std::string, std::string>> remOverrepBulk(std::vector<SRA> & sras, std::string threads, std::string ram_gb,
                    bool dispOutput, bool retainInterFiles, bool compressFiles,
                    std::string logFilePath, std::string outDir, bool noKrakenDbs,
                    const INI_MAP & cfgIni);
