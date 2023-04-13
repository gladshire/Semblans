#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <thread>
#include <boost/filesystem.hpp>
#include "sra.h"
#include "ini_parse.h"

namespace fs = boost::filesystem;


class transcript {
  private:
    std::string org_name;
    std::string tax_id;
    fs::path trans_path_trinity;
    fs::path trans_path_gene_map;
    fs::path trans_path_blastx;
    fs::path trans_path_ccut;
    fs::path trans_path_cinfo;
    fs::path trans_path_chimera;
    fs::path trans_path_index;
    fs::path trans_path_quant;
    fs::path trans_path_clust;
    fs::path trans_path_largest;
    fs::path trans_path_redund;
    fs::path trans_path_counts;
    fs::path trans_path_cds;
    fs::path trans_path_prot;
  public:
    transcript();
    transcript(SRA sra);
    transcript(std::string filename, INI_MAP cfgIni);
    std::string get_org_name();
    std::string get_tax_id();
    fs::path get_trans_path_trinity();
    fs::path get_trans_path_gene_map();
    fs::path get_trans_path_blastx();
    fs::path get_trans_path_ccut();
    fs::path get_trans_path_cinfo();
    fs::path get_trans_path_chimera();
    fs::path get_trans_path_index();
    fs::path get_trans_path_quant();
    fs::path get_trans_path_clust();
    fs::path get_trans_path_largest();
    fs::path get_trans_path_redund();
    fs::path get_trans_path_counts();
    fs::path get_trans_path_cds();
    fs::path get_trans_path_prot();
    void set_trans_path_trinity(fs::path transTrinity);
    void set_trans_path_blastx(fs::path transBlastx);
    void set_trans_path_ccut(fs::path transCcut);
    void set_trans_path_cinfo(fs::path transCinfo);
    void set_trans_path_chimera(fs::path transChimera);
    void set_trans_path_index(fs::path transIndex);
    void set_trans_path_quant(fs::path transQuant);
    void set_trans_path_clust(fs::path transClust);
    void set_trans_path_largest(fs::path transLargest);
    void set_trans_path_redund(fs::path transRedund);
    void set_trans_path_counts(fs::path transCounts);
    void set_trans_path_cds(fs::path transCds);
    void set_trans_path_prot(fs::path transProt);
    std::string make_file_str();
    std::string makeCheckpointName(std::string ext);
    void makeCheckpoint(std::string ext);
    bool checkpointExists(std::string ext);
};
