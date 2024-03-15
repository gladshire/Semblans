#pragma once
#include <iostream>
#include <iomanip>
#include <thread>
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
#include "annotate.h"
#include "print_info.h"


std::vector<transcript> get_transcript(std::vector<SRA> sras);

bool stringToBool(std::string boolStr);

void blastxBulk(const std::vector<transcript> & transVec, std::string threads,
                bool dispOutput, std::string logFilePath,
                std::string refProt, std::string outDir,
                const INI_MAP & cfgIni);

std::vector<std::string> remChimeraBulk(const std::vector<transcript> & transVec, std::string ram_gb,
                                        bool retainInterFiles, bool dispOutput, std::string logFilePath,
                                        const INI_MAP & cfgIni);

void salmonBulk(const std::vector<transcript> & transVec, std::string threads,
                bool retainInterFiles, bool dispOutput, std::string logFilePath,
                std::vector<std::string> readFiles1,
                std::vector<std::string> readFiles2,
                std::string outDir,
                const INI_MAP & cfgIni);

std::vector<std::string> corsetBulk(const std::vector<transcript> & transVec, std::string ram_gb,
                                    bool retainInterFiles, bool dispOutput, std::string logFilePath,
                                    const INI_MAP & cfgIni);

std::vector<std::pair<std::string, std::string>> transdecBulk(const std::vector<transcript> & transVec,
                                                              std::string threads, std::string ram_gb,
                                                              bool retainInterFiles, bool dispOutput,
                                                              std::string logFilePath,
                                                              std::string refProt, std::string outDir,
                                                              const INI_MAP & cfgIni);
