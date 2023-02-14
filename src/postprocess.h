#pragma once
#include <iostream>
#include <iomanip>
#include <sys/ioctl.h>
#include <unistd.h>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include <confini.h>
#include "sra.h"
#include "sra_toolkit.h"
#include "transcript.h"
#include "seq.h"
#include "ini_parse.h"
#include "ncbi_blast.h"
#include "rem_chimera.h"
#include "salmon_wrap.h"
#include "corset_wrap.h"
#include "filter_corset.h"
#include "transdecoder_wrap.h"
#include "print_info.h"

std::vector<transcript> get_transcript(std::vector<SRA> sras);

void print_help();
