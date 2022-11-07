#include "transcript.h"

/*
std::vector<std::string> stepDirs = {"00-Raw_reads/", "01-Quality_analysis_1/",
                                     "02-Error_correction/", "03-Trimming/",
                                     "04-Filter_foreign/", "05-Quality_analysis_2/",
                                     "06-Filter_overrepresented/", "07-Trinity_assembly/",
                                     "08-Filter_chimera/"};
*/
extern std::vector<std::string> stepDirs;

transcript::transcript() {
  org_name = "";
  tax_id = "";
}

transcript::transcript(SRA sra) {
  std::string projPath(sra.get_sra_path_raw().first.parent_path().parent_path().c_str());
  projPath += "/";
  org_name = sra.get_org_name();
  tax_id = sra.get_tax_id();
  std::string fileBase = make_file_str();
  
  trans_path_trinity = (projPath + stepDirs[7] + fileBase + ".Trinity.fasta").c_str();
  trans_path_gene_map = (projPath + stepDirs[7] + fileBase + "Trinity.fasta" + 
                         ".gene_trans_map").c_str();
  trans_path_chimera = (projPath + stepDirs[9] + fileBase + ".chim_filt.fasta").c_str();
  // Define salmon index path
  // Define salmon quant path
  // Define corset cluster path
  // Define transdecoder cds path
  // Define transdecoder prot path 
}

std::string transcript::get_org_name() {
  return org_name;
}

std::string transcript::get_tax_id() {
  return tax_id;
}

fs::path transcript::get_trans_path_trinity() {
  return trans_path_trinity;
}

fs::path transcript::get_trans_path_gene_map() {
  return trans_path_gene_map;
}

fs::path transcript::get_trans_path_chimera() {
  return trans_path_chimera;
}

fs::path transcript::get_trans_path_index() {
  return trans_path_index;
}

fs::path transcript::get_trans_path_quant() {
  return trans_path_quant;
}

fs::path transcript::get_trans_path_clust() {
  return trans_path_clust;
}

fs::path transcript::get_trans_path_cds() {
  return trans_path_cds;
}

fs::path transcript::get_trans_path_prot() {
  return trans_path_prot;
}

std::string transcript::make_file_str() {
  std::string filename;
  filename = get_tax_id() + "_" + get_org_name().replace(get_org_name().find(" "), 1, "_");
  return filename;
}
