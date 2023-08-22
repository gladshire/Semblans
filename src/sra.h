#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <chrono>
#include <thread>
#include <boost/filesystem.hpp>
#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_utils.hpp>
#include "ini_parse.h"
#include "log.h"

namespace fs = boost::filesystem;

// SRA class definition
// SRA run information is collated into single object
class SRA {
  private:
    std::string sra_accession;
    std::string org_name;
    std::string tax_id;
    long int spots;
    int spots_m;
    uintmax_t bp;
    bool paired;
    bool compressedFiles;
    std::string file_prefix_1;
    std::string file_prefix_2;
    fs::path fastqc_dir_1_1;
    fs::path fastqc_dir_2_1;
    fs::path sra_path_raw_1;
    fs::path sra_path_raw_2;
    fs::path sra_path_corr_1;
    fs::path sra_path_corr_2;
    fs::path sra_path_corr_fix_1;
    fs::path sra_path_corr_fix_2;
    fs::path sra_path_trim_u1;
    fs::path sra_path_trim_u2;
    fs::path sra_path_trim_p1;
    fs::path sra_path_trim_p2;
    fs::path sra_path_for_filt_1;
    fs::path sra_path_for_filt_2;
    fs::path fastqc_dir_1_2;
    fs::path fastqc_dir_2_2;
    fs::path sra_path_orep_filt_1;
    fs::path sra_path_orep_filt_2;
  public:
    SRA();
    SRA(std::string sra_accession, INI_MAP cfgIni, bool dispOutput,
        bool compressedFiles, std::string logFile, int num);
    SRA(std::string fileName1, std::string fileName2, INI_MAP cfgIni, bool compressedFiles);
    SRA(const SRA & sra);
    std::string get_accession() const;
    std::string get_org_name() const;
    std::string get_tax_id() const;
    long unsigned int get_num_reads() const;
    void set_num_reads(long unsigned int numReads);
    int get_spots_m() const;
    uintmax_t get_bp() const;
    bool is_paired() const;
    std::pair<std::string, std::string> get_file_prefix() const;
    std::pair<fs::path, fs::path> get_fastqc_dir_1() const;
    std::pair<fs::path, fs::path> get_sra_path_raw() const;
    std::pair<fs::path, fs::path> get_sra_path_corr() const;
    std::pair<fs::path, fs::path> get_sra_path_corr_fix() const;
    std::pair<fs::path, fs::path> get_sra_path_trim_u() const;
    std::pair<fs::path, fs::path> get_sra_path_trim_p() const;
    std::pair<fs::path, fs::path> get_sra_path_for_filt() const;
    std::pair<fs::path, fs::path> get_fastqc_dir_2() const;
    std::pair<fs::path, fs::path> get_sra_path_orep_filt() const;
    void set_num_reads(long int numReads);
    void set_fastqc_dir_1(std::pair<fs::path, fs::path> fastqcFiles1);
    void set_sra_path_raw(std::pair<fs::path, fs::path> sraRawFiles);
    void set_sra_path_corr(std::pair<fs::path, fs::path> sraCorrFiles);
    void set_sra_path_corr_fix(std::pair<fs::path, fs::path> sraCorrFixFiles);
    void set_sra_path_trim_u(std::pair<fs::path, fs::path> sraTrimFilesU);
    void set_sra_path_trim_p(std::pair<fs::path, fs::path> sraTrimFilesP);
    void set_sra_path_for_filt(std::pair<fs::path, fs::path> sraForFiltFiles);
    void set_fastqc_dir_2(std::pair<fs::path, fs::path> fastqcFiles2);
    void set_sra_path_orep_filt(std::pair<fs::path, fs::path> sraOrepFiltFiles);
    std::string make_file_str();
    std::string makeCheckpointName(std::string ext);
    void makeCheckpoint(std::string ext);
    bool checkpointExists(std::string ext);
};
