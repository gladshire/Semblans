#pragma once
#include <iostream>
#include <iomanip>
#include <sys/ioctl.h>
#include <unistd.h>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include "sra.h"
#include "sra_toolkit.h"
#include "transcript.h"
#include "seq.h"
#include "ini_parse.h"
#include "ncbi_blast.h"
#include "diamond.h"
#include "rem_chimera.h"
#include "salmon_wrap.h"
#include "corset_wrap.h"
#include "filter_corset.h"
#include "transdecoder_wrap.h"
#include "print_info.h"

std::vector<transcript> get_transcript(std::vector<SRA> sras);

void blastxDiamBulk(const std::vector<transcript> & transVec, std::string threads,
                    bool dispOutput, std::string logFilePath, const INI_MAP & cfgIni);

void remChimeraBulk(const std::vector<transcript> & transVec, std::string ram_gb,
                    std::string logFilePath);

void salmonBulk(const std::vector<transcript> & transVec, std::string threads,
                bool dispOutput, std::string logFilePath);

void corsetBulk(const std::vector<transcript> & transVec, std::string ram_gb,
                bool dispOutput, std::string logFilePath);

void transdecBulk(const std::vector<transcript> & transVec, std::string threads,
                  std::string ram_gb, bool dispOutput, std::string logFilePath,
                  const INI_MAP & cfgIni);
