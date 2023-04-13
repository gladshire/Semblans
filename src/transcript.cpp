#include "transcript.h"

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
 
  // Define trinity assembly path
  trans_path_trinity = (projPath + stepDirs[7] + fileBase + ".Trinity.fasta").c_str();
  // Define trinity gene map path
  trans_path_gene_map = (projPath + stepDirs[7] + fileBase + ".Trinity.fasta" + 
                         ".gene_trans_map").c_str();
  // Define chimera cut file path
  trans_path_ccut = (projPath + stepDirs[8] + fileBase + ".Trinity.cut").c_str();
  // Define chimera info file path
  trans_path_cinfo = (projPath + stepDirs[8] + fileBase + ".Trinity.info").c_str();
  // Define chimera-filtered transcript path
  trans_path_chimera = (projPath + stepDirs[8] + fileBase + ".chim_filt.fasta").c_str();
  // Define blastx output path
  trans_path_blastx = (projPath + stepDirs[8] + fileBase + ".Trinity.blastx").c_str();
  // Define salmon index path
  trans_path_index = (projPath + stepDirs[9] + fileBase + "_salmon_index").c_str();
  // Define salmon quant path
  trans_path_quant = (projPath + stepDirs[9] + fileBase + "_salmon_quant").c_str();
  // Define corset cluster path
  trans_path_clust = (projPath + stepDirs[9] + fileBase + "_salmon-clusters.txt").c_str();
  // Define largest cluster transcripts path
  trans_path_largest = (projPath + stepDirs[9] + fileBase + ".lc.fasta").c_str();
  // Define redundant transcripts path
  trans_path_redund = (projPath + stepDirs[9] + fileBase + ".rd.fasta").c_str();
  // Define corset counts path
  trans_path_counts = (projPath + stepDirs[9] + fileBase + "_salmon-counts.txt").c_str();
  // Define transdecoder cds path
  trans_path_cds = (projPath + stepDirs[10] + ".transdecoder.cds").c_str();
  // Define transdecoder prot path 
  trans_path_prot = (projPath + stepDirs[10] + ".transdecoder.pep").c_str();

} 

transcript::transcript(std::string filename, INI_MAP cfgIni) {
  std::string projPath(fs::path(fs::path(cfgIni["General"]["output_directory"].c_str()) /
                                fs::path(cfgIni["General"]["project_name"].c_str())).c_str());
  org_name = "";
  tax_id = "";
  std::string fileBase(fs::path(filename).filename().stem().c_str());
  
  // Define trinity assembly path
  trans_path_trinity = (projPath + stepDirs[7] + fileBase + ".Trinity.fasta").c_str();
  // Define trinity gene map path
  trans_path_gene_map = (projPath + stepDirs[7] + fileBase + ".Trinity.fasta" + 
                         ".gene_trans_map").c_str();
  // Define chimera cut file path
  trans_path_ccut = (projPath + stepDirs[8] + fileBase + ".Trinity.cut").c_str();
  // Define chimera info file path
  trans_path_cinfo = (projPath + stepDirs[8] + fileBase + ".Trinity.info").c_str();
  // Define chimera-filtered transcript path
  trans_path_chimera = (projPath + stepDirs[8] + fileBase + ".chim_filt.fasta").c_str();
  // Define blastx output path
  trans_path_blastx = (projPath + stepDirs[8] + fileBase + ".Trinity.blastx").c_str();
  // Define salmon index path
  trans_path_index = (projPath + stepDirs[9] + fileBase + "_salmon_index").c_str();
  // Define salmon quant path
  trans_path_quant = (projPath + stepDirs[9] + fileBase + "_salmon_quant").c_str();
  // Define corset cluster path
  trans_path_clust = (projPath + stepDirs[9] + fileBase + "_salmon-clusters.txt").c_str();
  // Define largest cluster transcripts path
  trans_path_largest = (projPath + stepDirs[9] + fileBase + ".lc.fasta").c_str();
  // Define redundant transcripts path
  trans_path_redund = (projPath + stepDirs[9] + fileBase + ".rd.fasta").c_str();
  // Define corset counts path
  trans_path_counts = (projPath + stepDirs[9] + fileBase + "_salmon-counts.txt").c_str();
  // Define transdecoder cds path
  trans_path_cds = (projPath + stepDirs[10] + ".transdecoder.cds").c_str();
  // Define transdecoder prot path 
  trans_path_prot = (projPath + stepDirs[10] + ".transdecoder.pep").c_str();
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

fs::path transcript::get_trans_path_ccut() {
  return trans_path_ccut;
}

fs::path transcript::get_trans_path_cinfo() {
  return trans_path_cinfo;
}

fs::path transcript::get_trans_path_chimera() {
  return trans_path_chimera;
}

fs::path transcript::get_trans_path_blastx() {
  return trans_path_blastx;
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

fs::path transcript::get_trans_path_largest() {
  return trans_path_largest;
}

fs::path transcript::get_trans_path_redund() {
  return trans_path_redund;
}

fs::path transcript::get_trans_path_counts() {
  return trans_path_counts;
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

std::string transcript::makeCheckpointName(std::string ext) {
  fs::path outDir = get_trans_path_trinity().parent_path().parent_path() / "checkpoints";
  std::string cpFileName = std::string(outDir.c_str()) + "/" + make_file_str() + "." + ext + ".ok";
  return cpFileName;
}

void transcript::makeCheckpoint(std::string ext) {
  std::string cpFileName = makeCheckpointName(ext);
  std::ofstream cpFile;
  cpFile.open(cpFileName);
  cpFile.close();
}

bool transcript::checkpointExists(std::string ext) {
  fs::path cpFilePath(makeCheckpointName(ext).c_str());
  if (fs::exists(cpFilePath)) {
    return true;
  }
  else {
    return false;
  }
}
