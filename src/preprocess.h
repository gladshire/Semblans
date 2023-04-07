#pragma once
#include <iostream>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sys/ioctl.h>
#include <unistd.h>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
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

void retrieve_sra_data(const std::vector<SRA> & sras, std::string threads,
                       bool dispOutput, std::string logFile);

bool stringToBool(std::string boolStr);

void print_help();

void fastqcBulk1(const std::vector<SRA> & sras, std::string threads, bool dispOutput,
                 std::string logFilePath);

void fastqcBulk2(const std::vector<SRA> & sras, std::string threads, bool dispOutput,
                 std::string logFilePath);

void errorCorrBulk(const std::vector<SRA> & sras, std::string threads, bool dispOutput,
                   std::string logFilePath);

void remUnfixBulk(const std::vector<SRA> & sras, std::string threads, std::string ram_gb,
                  bool dispOutput, std::string logFilePath);

void trimBulk(const std::vector<SRA> & sras, std::string threads, bool dispOutput,
              std::string logFilePath);

void filtForeignBulk(const std::vector<SRA> & sras, std::vector<std::string> krakenDbs,
                     std::string krakenConf, std::string threads, bool dispOutput,
                     std::string logFilePath);

void remOverrepBulk(const std::vector<SRA> & sras, std::string threads, std::string ram_gb,
                    bool dispOutput, std::string logFilePath);
