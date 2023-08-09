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
#include "print_info.h"

#define PATH_PIGZ std::string((dl::program_location().parent_path() / fs::path(std::string("../external/pigz/pigz"))).c_str())


void progressAnim(int numSpace);

void retrieve_sra_data(const std::vector<SRA> & sras, std::string threads,
                       bool dispOutput, bool retainInterFiles,
                       bool compressOutput, std::string logFile);

bool stringToBool(std::string boolStr);

void print_help();

void fastqcBulk1(const std::vector<SRA> & sras, std::string threads, bool dispOutput,
                 std::string logFilePath);

void fastqcBulk2(const std::vector<SRA> & sras, std::string threads, bool dispOutput,
                 std::string logFilePath, const INI_MAP & cfgIni);

void errorCorrBulk(const std::vector<SRA> & sras, std::string threads,
                   bool dispOutput, bool retainInterFiles, bool compressFiles,
                   std::string logFilePath, const INI_MAP & cfgIni);

bool remUnfixBulk(const std::vector<SRA> & sras, std::string threads, std::string ram_gb,
                  bool dispOutput, bool retainInterFiles, bool compressFiles,
                  std::string logFilePath, const INI_MAP & cfgIni);

void trimBulk(const std::vector<SRA> & sras, std::string threads,
              bool dispOutput, bool retainInterFiles,
              std::string logFilePath, const INI_MAP & cfgIni);

void filtForeignBulk(const std::vector<SRA> & sras, std::vector<std::string> krakenDbs,
                     std::string threads, bool dispOutput, bool compressFiles, bool retainInterFiles,
                     std::string logFilePath, const INI_MAP & cfgIni);

bool remOverrepBulk(const std::vector<SRA> & sras, std::string threads, std::string ram_gb,
                    bool dispOutput, bool retainInterFiles, bool compressFiles,
                    std::string logFilePath, const INI_MAP & cfgIni);
