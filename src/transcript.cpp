#include "transcript.h"

extern std::vector<std::string> stepDirs;

// Default constructor for transcript object
transcript::transcript() {
  org_name = "";
  tax_id = "";
}

// Constructor for transcript object, taking an SRA run object as input
transcript::transcript(SRA sra) {
  //std::string projPath(sra.get_sra_path_trim_u().first.parent_path().parent_path().c_str());
  std::string projPath(sra.get_fastqc_dir_2().first.parent_path().parent_path().parent_path().c_str());
  projPath += "/";
  org_name = sra.get_org_name();
  tax_id = sra.get_tax_id();
  std::string fileBase;
  if (sra.get_accession() != "") {
    fileBase = sra.make_file_str();
  }
  else {
    fileBase = sra.get_file_prefix().first;
  }
  file_prefix = fileBase;
 
  // Define trinity assembly path
  trans_path_trinity = (projPath + stepDirs[7] + "/" + fileBase + ".Trinity.fasta").c_str();
  // Define trinity gene map path
  trans_path_gene_map = (projPath + stepDirs[7] + "/" + fileBase + ".Trinity.fasta" + 
                         ".gene_trans_map").c_str();
  // Define chimera cut file path
  trans_path_ccut = (projPath + stepDirs[8] + "/" + fileBase + ".Trinity.cut").c_str();
  // Define chimera info file path
  trans_path_cinfo = (projPath + stepDirs[8] + "/" + fileBase + ".Trinity.info").c_str();
  // Define chimera-filtered transcript path
  trans_path_chimera = (projPath + stepDirs[8] + "/" + fileBase + ".chim_filt.fasta").c_str();
  // Define blastx output path
  trans_path_blastx = (projPath + stepDirs[8] + "/" + fileBase + ".Trinity.blastx").c_str();
  // Define salmon index path
  trans_path_index = (projPath + stepDirs[9] + "/" + fileBase + "_salmon_index").c_str();
  // Define salmon quant path
  trans_path_quant = (projPath + stepDirs[9] + "/" + fileBase + "_salmon_quant").c_str();
  // Define corset cluster path
  trans_path_clust = (projPath + stepDirs[9] + "/" + fileBase + "_salmon-clusters.txt").c_str();
  // Define largest cluster transcripts path
  trans_path_largest = (projPath + stepDirs[9] + "/" + fileBase + ".lc.fasta").c_str();
  // Define redundant transcripts path
  trans_path_redund = (projPath + stepDirs[9] + "/" + fileBase + ".rd.fasta").c_str();
  // Define corset counts path
  trans_path_counts = (projPath + stepDirs[9] + "/" + fileBase + "_salmon-counts.txt").c_str();
  // Define transdecoder cds path
  trans_path_cds = (projPath + stepDirs[10] + "/" + fileBase + ".transdecoder.cds.fasta").c_str();
  // Define transdecoder prot path 
  trans_path_prot = (projPath + stepDirs[10] + "/" + fileBase + ".transdecoder.pep.fasta").c_str();
  // Define annotated transcripts path
  trans_path_annot = (projPath + stepDirs[11] + "/" + fileBase + ".annotated.fasta").c_str();
} 

// Constructor for transcript object, taking a FASTA file as input
transcript::transcript(std::string filename, INI_MAP cfgIni) {
  std::string projPath(fs::canonical(fs::path(cfgIni["General"]["output_directory"].c_str()) /
                                     fs::path(cfgIni["General"]["project_name"].c_str())).c_str());
  projPath += "/";
  org_name = "";
  tax_id = "";
  std::string fileBase(fs::path(filename).filename().stem().stem().c_str());
  //while (fileBase.find(".") != std::string::npos) {
  //  fileBase = std::string(fs::path(fileBase).stem().c_str());
  //}
  file_prefix = fileBase;
  
  // Define trinity assembly path
  trans_path_trinity = (projPath + stepDirs[7] + "/" + fileBase + ".Trinity.fasta").c_str();
  // Define trinity gene map path
  trans_path_gene_map = (projPath + stepDirs[7] + "/" + fileBase + ".Trinity.fasta" + 
                         ".gene_trans_map").c_str();
  // Define chimera cut file path
  trans_path_ccut = (projPath + stepDirs[8] + "/" + fileBase + ".Trinity.cut").c_str();
  // Define chimera info file path
  trans_path_cinfo = (projPath + stepDirs[8] + "/" + fileBase + ".Trinity.info").c_str();
  // Define chimera-filtered transcript path
  trans_path_chimera = (projPath + stepDirs[8] + "/" + fileBase + ".chim_filt.fasta").c_str();
  // Define blastx output path
  trans_path_blastx = (projPath + stepDirs[8] + "/" + fileBase + ".Trinity.blastx").c_str();
  // Define salmon index path
  trans_path_index = (projPath + stepDirs[9] + "/" + fileBase + "_salmon_index").c_str();
  // Define salmon quant path
  trans_path_quant = (projPath + stepDirs[9] + "/" + fileBase + "_salmon_quant").c_str();
  // Define corset cluster path
  trans_path_clust = (projPath + stepDirs[9] + "/" + fileBase + "_salmon-clusters.txt").c_str();
  // Define largest cluster transcripts path
  trans_path_largest = (projPath + stepDirs[9] + "/" + fileBase + ".lc.fasta").c_str();
  // Define redundant transcripts path
  trans_path_redund = (projPath + stepDirs[9] + "/" + fileBase + ".rd.fasta").c_str();
  // Define corset counts path
  trans_path_counts = (projPath + stepDirs[9] + "/" + fileBase + "_salmon-counts.txt").c_str();
  // Define transdecoder cds path
  trans_path_cds = (projPath + stepDirs[10] + "/" + fileBase + ".transdecoder.cds.fasta").c_str();
  // Define transdecoder prot path 
  trans_path_prot = (projPath + stepDirs[10] + "/" + fileBase + ".transdecoder.pep.fasta").c_str();
  // Define annotated transcripts path
  trans_path_annot = (projPath + stepDirs[11] + "/" + fileBase + ".annotated.fasta").c_str();
}

transcript::transcript(std::string filename, std::string outDir) {
  outDir += "/";
  org_name = "";
  tax_id = "";
  fs::path filePrefix(fs::path(filename.c_str()));
  while (!filePrefix.extension().empty()) {
    filePrefix = filePrefix.stem();
  }
  //while (fileBase.find(".") != std::string::npos) {
  //  fileBase = std::string(fs::path(fileBase).stem().c_str());
  //}
  std::string fileBase = std::string(filePrefix.c_str());

  file_prefix = fileBase;
  
  // Define trinity assembly path
  //trans_path_trinity = (projPath + stepDirs[7] + "/" + fileBase + ".Trinity.fasta").c_str();
  trans_path_trinity = filename;
  // Define trinity gene map path
  //trans_path_gene_map = (projPath + stepDirs[7] + "/" + fileBase + ".Trinity.fasta" + 
  //                       ".gene_trans_map").c_str();
  // Define chimera cut file path
  trans_path_ccut = (outDir + stepDirs[8] + "/" + fileBase + ".cut").c_str();
  // Define chimera info file path
  trans_path_cinfo = (outDir + stepDirs[8] + "/" + fileBase + ".info").c_str();
  // Define chimera-filtered transcript path
  trans_path_chimera = (outDir + stepDirs[8] + "/" + fileBase + ".chim_filt.fasta").c_str();
  // Define blastx output path
  trans_path_blastx = (outDir + stepDirs[8] + "/" + fileBase + ".blastx").c_str();
  // Define salmon index path
  trans_path_index = (outDir + stepDirs[9] + "/" + fileBase + "_salmon_index").c_str();
  // Define salmon quant path
  trans_path_quant = (outDir + stepDirs[9] + "/" + fileBase + "_salmon_quant").c_str();
  // Define corset cluster path
  trans_path_clust = (outDir + stepDirs[9] + "/" + fileBase + "_salmon-clusters.txt").c_str();
  // Define largest cluster transcripts path
  trans_path_largest = (outDir + stepDirs[9] + "/" + fileBase + ".lc.fasta").c_str();
  // Define redundant transcripts path
  trans_path_redund = (outDir + stepDirs[9] + "/" + fileBase + ".rd.fasta").c_str();
  // Define corset counts path
  trans_path_counts = (outDir + stepDirs[9] + "/" + fileBase + "_salmon-counts.txt").c_str();
  // Define transdecoder cds path
  trans_path_cds = (outDir + stepDirs[10] + "/" + fileBase + ".transdecoder.cds.fasta").c_str();
  // Define transdecoder prot path 
  trans_path_prot = (outDir + stepDirs[10] + "/" + fileBase + ".transdecoder.pep.fasta").c_str();
  // Define annotated transcripts path
  trans_path_annot = (outDir + stepDirs[11] + "/" + fileBase + ".annotated.fasta").c_str();
}

// Return the organism name
std::string transcript::get_org_name() {
  return org_name;
}

// Return the taxon ID
std::string transcript::get_tax_id() {
  return tax_id;
}

// Return the path to the Trinity assembly file
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

fs::path transcript::get_trans_path_annot() {
  return trans_path_annot;
}

std::string transcript::get_file_prefix() {
  return file_prefix;
}

void transcript::set_trans_path_trinity(fs::path transTrinity) {
  trans_path_trinity = transTrinity;
}

void transcript::set_trans_path_blastx(fs::path transBlastx) {
  trans_path_blastx = transBlastx;
}

void transcript::set_trans_path_ccut(fs::path transCcut) {
  trans_path_ccut = transCcut;
}

void transcript::set_trans_path_cinfo(fs::path transCinfo) {
  trans_path_cinfo = transCinfo;
}

void transcript::set_trans_path_chimera(fs::path transChimera) {
  trans_path_chimera = transChimera;
}

void transcript::set_trans_path_index(fs::path transIndex) {
  trans_path_index = transIndex;
}

void transcript::set_trans_path_quant(fs::path transQuant) {
  trans_path_quant = transQuant;
}

void transcript::set_trans_path_clust(fs::path transClust) {
  trans_path_clust = transClust;
}

void transcript::set_trans_path_largest(fs::path transLargest) {
  trans_path_largest = transLargest;
}

void transcript::set_trans_path_redund(fs::path transRedund) {
  trans_path_redund = transRedund;
}

void transcript::set_trans_path_counts(fs::path transCounts) {
  trans_path_counts = transCounts;
}

void transcript::set_trans_path_cds(fs::path transCds) {
  trans_path_cds = transCds;
}

void transcript::set_trans_path_prot(fs::path transProt) {
  trans_path_prot = transProt;
}

std::string transcript::make_file_str() {
  std::string filename;
  filename = get_tax_id() + "_" + get_org_name().replace(get_org_name().find(" "), 1, "_");
  return filename;
}

std::string transcript::makeCheckpointName(std::string ext) {
  fs::path outDir = get_trans_path_trinity().parent_path().parent_path() / ".checkpoints";
  std::string cpFileName = std::string(outDir.c_str()) + "/" + get_file_prefix() + "." + ext + ".ok";
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
