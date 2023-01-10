#pragma once
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/dll.hpp>
#include <confini.h>
#include "sra.h"
#include "sra_toolkit.h"
#include "transcript.h"
#include "ini_parse.h"
#include "trinity_wrap.h"

std::vector<transcript> get_transcript(std::vector<SRA> sras);

void print_help();
