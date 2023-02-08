#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include "sra.h"
#include "transcript.h"
#include "seq.h"
#include "seq_hash.h"
#include "print_info.h"

namespace fs = boost::filesystem;
namespace dl = boost::filesystem;


void filterCorset(transcript trans, std::string clusterPath, std::string out_dir);
