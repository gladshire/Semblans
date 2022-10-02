#include "sra.h"
#include "ini_parse.h"


// Default constructor for SRA object
SRA::SRA() {
  sra_accession = "";
  org_name = "";
  tax_id = "";
  spots = -1;
  spots_m = -1;
  bp = -1;
  paired = false;
}

// SRA accession number constructor for SRA object
// Takes SRA accession number string as input
// Fills object members with correct values using NCBI eutils API
SRA::SRA(std::string sra_accession) {
  extern INI_MAP cfgIni;
  std::string outDir(cfgIni["General"]["output_directory"]);
  std::string projName(cfgIni["General"]["project_name"]);

  // Download temp XML file for SRA accession, containing information for object members
  std::string curlCmdStr = "curl -s -o tmp.xml \"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?rettype=runinfo&db=sra&id=\"" + sra_accession;
  system(curlCmdStr.c_str());
  
  // Parse XML file for object information
  rapidxml::file<> sra_xml("tmp.xml");
  rapidxml::xml_document<> sra_doc;
  sra_doc.parse<0>(sra_xml.data());

  rapidxml::xml_node<> * parse_node = sra_doc.first_node()->first_node()->first_node();
 
  // Set object SRA accession number
  this->sra_accession = sra_accession;

  // Set object member for number of spots read
  while(strcmp(parse_node->name(), "spots")) {
    parse_node = parse_node->next_sibling();
  }
  spots = atoi(parse_node->value());

  // Set object member for number of base pairs in SRA run
  while(strcmp(parse_node->name(), "bases")) {
    parse_node = parse_node->next_sibling();
  }
  bp = atoi(parse_node->value());

  // Set object member for number of spots that have mates
  while(strcmp(parse_node->name(), "spots_with_mates")) {
    parse_node = parse_node->next_sibling();
  }
  spots_m = atoi(parse_node->value());

  // Set object member for whether SRA run is single-ended or paired-ended
  while(strcmp(parse_node->name(), "LibraryLayout")) {
    parse_node = parse_node->next_sibling();
  }
  if (!strcmp(parse_node->value(), "PAIRED")) {
    paired = true;
  }
  else {
    paired = false;
  }

  // Set object member for taxonomic ID number
  while(strcmp(parse_node->name(), "TaxID")) {
    parse_node = parse_node->next_sibling();
  }
  this->tax_id = std::string(parse_node->value());

  // Set object member for organism name
  while(strcmp(parse_node->name(), "ScientificName")) {
    parse_node = parse_node->next_sibling();
  }
  this->org_name = std::string(parse_node->value());

  std::string fileBase = make_file_str();
  std::string projPath = outDir + projName + "/";
  std::vector<std::string> stepDirs = {"00-Raw_Reads/", "01-Quality_analysis_1/",
                                       "02-Error_correction/", "03-Trimming/",
                                       "04-Filtering/"};
  sra_path_raw_1  = (projPath + stepDirs[0] + fileBase + ".fastq").c_str();
  fastqc_dir_1    = (projPath + stepDirs[1] + fileBase + "/" + fileBase).c_str();
  sra_path_corr_1 = (projPath + stepDirs[2] + fileBase + ".corr.fastq").c_str();
  sra_path_trim_1 = (projPath + stepDirs[3] + fileBase + ".trim.fastq").c_str();
  sra_path_filt_1 = (projPath + stepDirs[4] + fileBase + ".filt.fastq").c_str();
 
  if (paired) {
    std::string sra_path_raw_1_str(sra_path_raw_1.native().c_str());
    std::string fastqc_dir_1_str(fastqc_dir_1.native().c_str());
    std::string sra_path_corr_1_str(sra_path_corr_1.native().c_str());
    std::string sra_path_trim_1_str(sra_path_trim_1.native().c_str());
    std::string sra_path_filt_1_str(sra_path_filt_1.native().c_str());
    sra_path_raw_1 = (sra_path_raw_1_str.insert(sra_path_raw_1_str.length() - 6, "_1")).c_str();
    fastqc_dir_1 = (fastqc_dir_1_str.insert(fastqc_dir_1_str.length() - 1, "_1")).c_str();
    sra_path_corr_1 = (sra_path_corr_1_str.insert(sra_path_corr_1_str.length() - 10, "_1")).c_str();
    sra_path_trim_1 = (sra_path_trim_1_str.insert(sra_path_trim_1_str.length() - 10, "_1")).c_str();
    sra_path_filt_1 = (sra_path_trim_1_str.insert(sra_path_filt_1_str.length() - 10, "_1")).c_str();

    sra_path_raw_2 = (projPath + stepDirs[0] + fileBase + "_2.fastq").c_str();
    fastqc_dir_2 = (projPath + stepDirs[1] + fileBase + "/" + fileBase + "_2").c_str();
    sra_path_corr_2 = (projPath + stepDirs[2] + fileBase + "_2.corr.fastq").c_str();
    sra_path_trim_2 = (projPath + stepDirs[3] + fileBase + "_2.trim.fastq").c_str();
    sra_path_filt_2 = (projPath + stepDirs[4] + fileBase + "_2.filt.fastq").c_str();
  }
 
  system("rm tmp.xml");
}

// Getter function for object SRA accession number
std::string SRA::get_accession() {
  return sra_accession;
}

// Getter function for object organism name
std::string SRA::get_org_name() {
  return org_name;
}

// Getter function for object taxonomic ID
std::string SRA::get_tax_id() {
  return tax_id;
}

// Getter function for object number of read spots
int SRA::get_spots() {
  return spots;
}

// Getter function for object number of read spots that are paired
int SRA::get_spots_m() {
  return spots_m;
}

// Getter function for object number of base pairs in SRA run
long int SRA::get_bp() {
  return bp;
}

// Getter function for object single-end or paired-end layout
bool SRA::is_paired() {
  return paired;
}

// Getter function for path to SRA FastQC analysis directory
std::pair<fs::path, fs::path> SRA::get_fastqc_dir() {
  std::pair<fs::path, fs::path> fastqc_dir(fastqc_dir_1, fastqc_dir_2);
  return fastqc_dir;
}

// Getter function for path to raw SRA reads sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_raw() {
  std::pair<fs::path, fs::path> sra_path_raw(sra_path_raw_1, sra_path_raw_2);
  return sra_path_raw;
}

// Getter function for path to error-corrected SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_corr() {
  std::pair<fs::path, fs::path> sra_path_corr(sra_path_corr_1, sra_path_corr_2);
  return sra_path_corr;
}

// Getter function for path to trimmed SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_trim() {
  std::pair<fs::path, fs::path> sra_path_trim(sra_path_trim_1, sra_path_trim_2);
  return sra_path_trim;
}

// Getter function for path to filtered SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_filt() {
  std::pair<fs::path, fs::path> sra_path_filt(sra_path_filt_1, sra_path_filt_2);
  return sra_path_filt;
}

// Utility function to construct filename for SRA
std::string SRA::make_file_str() {
  std::string filename;
  filename = get_accession() + "_" +
             (std::string)(get_tax_id()) + "_" +
             get_org_name().replace(get_org_name().find(" "), 1, "_");
  return filename;
}

