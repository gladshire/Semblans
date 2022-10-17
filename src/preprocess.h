#pragma once
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/variant.hpp>
#include <confini.h>
#include "sra.h"
#include "ini_parse.h"
#include "sra_toolkit.h"
#include "fastqc_wrap.h"
#include "rcorr_wrap.h"
#include "rem_unfixable.h"
#include "trimm_wrap.h"
#include "kraken_wrap.h"
#include "print_info.h"

std::vector<SRA> get_sras(const INI_MAP &iniFile);

