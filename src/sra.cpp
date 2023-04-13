#include "sra.h"



std::vector<std::string> stepDirs = {"00-Raw_reads/", "01-Quality_analysis_1/",
                                     "02-Error_correction/", "03-Trimming/",
                                     "04-Filter_foreign/", "05-Quality_analysis_2/",
                                     "06-Filter_overrepresented/", "07-Trinity_assembly/",
                                     "08-Filter_chimera/", "09-Clustering/",
                                     "10-Final_cds_pep/"};


extern std::vector<std::string> stepDirs;

// Default constructor for SRA object
SRA::SRA() {
  sra_accession = "";
  org_name = "";
  tax_id = "";
  spots = -1;
  spots_m = -1;
  bp = -1;
  paired = false;
  file_prefix_1 = "";
  file_prefix_2 = "";
}

// SRA accession number constructor for SRA object
// Takes SRA accession number string as input
// Fills object members with correct values using NCBI eutils API
SRA::SRA(std::string sra_accession, INI_MAP cfgIni) {
  std::string outDir(cfgIni["General"]["output_directory"]);
  std::string projName(cfgIni["General"]["project_name"]);
  std::string apiKey(cfgIni["General"]["ncbi_api_key"]);
  std::chrono::milliseconds queryLim(500);
  if (apiKey != "") {
    std::chrono::milliseconds queryLim(250);
  }

  // Download temp XML file for SRA accession, containing information for object members
  std::string curlCmdStr = "curl -s -o tmp.xml \"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?rettype=runinfo&db=sra&id=\"" +
                           sra_accession + "&api_key=" + apiKey;
  int result;
  while (true) {
    system(curlCmdStr.c_str());
    std::this_thread::sleep_for(queryLim);
  
    // Parse XML file for object information
    try {
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
      break;
    }
    catch (const std::runtime_error & e) {
      std::cout << "Throttled: too many requests. Retrying ..." << std::endl;
      std::this_thread::sleep_for(5 * queryLim);
    }
  }
  
  std::string fileBase = make_file_str();
  file_prefix_1 = fileBase + "_1";
  file_prefix_2 = fileBase + "_2";

  std::string projPath = outDir + projName + "/";
  extern std::vector<std::string> stepDirs;

  sra_path_raw_1  = (projPath + stepDirs[0] + fileBase + ".fastq").c_str();
  fastqc_dir_1_1  = (projPath + stepDirs[1] + fileBase + "/" + fileBase).c_str();
  sra_path_corr_1 = (projPath + stepDirs[2] + fileBase + ".cor.fq").c_str();
  sra_path_corr_fix_1 = (projPath + stepDirs[2] + fileBase + ".cor.fix.fq").c_str();
  sra_path_trim_u1 = (projPath + stepDirs[3] + fileBase + ".trim.fq").c_str();
  sra_path_trim_p1 = (projPath + stepDirs[3] + fileBase + ".trim.fq").c_str();
  sra_path_for_filt_1 = (projPath + stepDirs[4] + fileBase + ".filt.fq").c_str();
  fastqc_dir_1_2 = (projPath + stepDirs[5] + fileBase + "/" + fileBase).c_str();
  sra_path_orep_filt_1 = (projPath + stepDirs[6] + fileBase + ".orep.filt.fq").c_str();
 
  if (paired) {
    std::string sra_path_raw_1_str(sra_path_raw_1.c_str());
    std::string fastqc_dir_1_1_str(fastqc_dir_1_1.c_str());
    std::string sra_path_corr_1_str(sra_path_corr_1.c_str());
    std::string sra_path_corr_fix_1_str(sra_path_corr_fix_1.c_str());
    std::string sra_path_trim_u1_str(sra_path_trim_u1.c_str());
    std::string sra_path_trim_p1_str(sra_path_trim_u1.c_str());
    std::string sra_path_for_filt_1_str(sra_path_for_filt_1.c_str());
    std::string fastqc_dir_1_2_str(fastqc_dir_1_2.c_str());
    std::string sra_path_orep_filt_1_str(sra_path_orep_filt_1.c_str());

    sra_path_raw_1 = (sra_path_raw_1_str.insert(sra_path_raw_1_str.length() - 6, "_1")).c_str();
    fastqc_dir_1_1 = (fastqc_dir_1_1_str.insert(fastqc_dir_1_1_str.length(), "_1")).c_str();
    sra_path_corr_1 = (sra_path_corr_1_str.insert(sra_path_corr_1_str.length() - 7, "_1")).c_str();
    sra_path_corr_fix_1 = (sra_path_corr_fix_1_str.insert(sra_path_corr_fix_1_str.length() - 11, "_1")).c_str();
    sra_path_trim_u1 = (sra_path_trim_u1_str.insert(sra_path_trim_u1_str.length() - 8, "_1.unpaired")).c_str();
    sra_path_trim_p1 = (sra_path_trim_p1_str.insert(sra_path_trim_p1_str.length() - 8, "_1.paired")).c_str();
    sra_path_for_filt_1 = (sra_path_for_filt_1_str.insert(sra_path_for_filt_1_str.length() - 8, "_1")).c_str();
    fastqc_dir_1_2 = (fastqc_dir_1_2_str.insert(fastqc_dir_1_2_str.length(), "_1")).c_str();
    sra_path_orep_filt_1 = (sra_path_orep_filt_1_str.insert(sra_path_orep_filt_1_str.length() - 13, "_1")).c_str();

    sra_path_raw_2 = (projPath + stepDirs[0] + fileBase + "_2.fastq").c_str();
    fastqc_dir_2_1 = (projPath + stepDirs[1] + fileBase + "/" + fileBase + "_2").c_str();
    sra_path_corr_2 = (projPath + stepDirs[2] + fileBase + "_2.cor.fq").c_str();
    sra_path_corr_fix_2 = (projPath + stepDirs[2] + fileBase + "_2.cor.fix.fq").c_str();
    sra_path_trim_u2 = (projPath + stepDirs[3] + fileBase + "_2.unpaired.trim.fq").c_str();
    sra_path_trim_p2 = (projPath + stepDirs[3] + fileBase + "_2.paired.trim.fq").c_str();
    sra_path_for_filt_2 = (projPath + stepDirs[4] + fileBase + "_2.filt.fq").c_str();
    fastqc_dir_2_2 = (projPath + stepDirs[5] + fileBase + "/" + fileBase + "_2").c_str();
    sra_path_orep_filt_2 = (projPath + stepDirs[6] + fileBase + "_2.orep.filt.fq").c_str();
  }
 
  system("rm tmp.xml");
}


// Constructor for if local paired-end data used
SRA::SRA(std::string fileName1, std::string fileName2, INI_MAP cfgIni) {
  std::string outDir(cfgIni["General"]["output_directory"]);
  std::string projName(cfgIni["General"]["project_name"]);
  std::string projPath = outDir + projName + "/";
  std::string localDataDir(cfgIni["General"]["local_data_directory"]);
  std::string fileBase1(fs::path(fileName1.c_str()).stem().c_str());
  std::string fileBase2(fs::path(fileName2.c_str()).stem().c_str());

  file_prefix_1 = fileBase1;
  file_prefix_2 = fileBase2;

  sra_accession = "";
  org_name = "";
  tax_id = "";
  spots = -1;
  spots_m = -1;
  bp = -1;
  paired = (fileName2 == "") ? false : true;
  
  //sra_path_raw_1  = (projPath + stepDirs[0] + fileBase1 + ".fastq").c_str();
  sra_path_raw_1 = (localDataDir + fileBase1 + ".fastq").c_str();
  fastqc_dir_1_1  = (projPath + stepDirs[1] + fileBase1 + "/" + fileBase1).c_str();
  sra_path_corr_1 = (projPath + stepDirs[2] + fileBase1 + ".cor.fq").c_str();
  sra_path_corr_fix_1 = (projPath + stepDirs[2] + fileBase1 + ".cor.fix.fq").c_str();
  sra_path_trim_u1 = (projPath + stepDirs[3] + fileBase1 + ".trim.fq").c_str();
  sra_path_trim_p1 = (projPath + stepDirs[3] + fileBase1 + ".trim.fq").c_str();
  sra_path_for_filt_1 = (projPath + stepDirs[4] + fileBase1 + ".filt.fq").c_str();
  fastqc_dir_1_2 = (projPath + stepDirs[5] + fileBase1 + "/" + fileBase1).c_str();
  sra_path_orep_filt_1 = (projPath + stepDirs[6] + fileBase1 + ".orep.filt.fq").c_str();

  if (paired) {
    //sra_path_raw_2 = (projPath + stepDirs[0] + fileBase2 + ".fastq").c_str();
    sra_path_raw_2 = (localDataDir + fileBase2 + ".fastq").c_str();
    fastqc_dir_2_1 = (projPath + stepDirs[1] + fileBase2 + "/" + fileBase2).c_str();
    sra_path_corr_2 = (projPath + stepDirs[2] + fileBase2 + ".cor.fq").c_str();
    sra_path_corr_fix_2 = (projPath + stepDirs[2] + fileBase2 + ".cor.fix.fq").c_str();
    sra_path_trim_u1 = (projPath + stepDirs[3] + fileBase1 + ".unpaired.trim.fq").c_str();
    sra_path_trim_p1 = (projPath + stepDirs[3] + fileBase1 + ".paired.trim.fq").c_str();
    sra_path_trim_u2 = (projPath + stepDirs[3] + fileBase2 + ".unpaired.trim.fq").c_str();
    sra_path_trim_p2 = (projPath + stepDirs[3] + fileBase2 + ".paired.trim.fq").c_str();
    sra_path_for_filt_2 = (projPath + stepDirs[4] + fileBase2 + ".filt.fq").c_str();
    fastqc_dir_2_2 = (projPath + stepDirs[5] + fileBase2 + "/" + fileBase2).c_str();
    sra_path_orep_filt_2 = (projPath + stepDirs[6] + fileBase2 + ".orep.filt.fq").c_str();
  }
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

std::pair<std::string, std::string> SRA::get_file_prefix() {
  std::pair<std::string, std::string> file_prefix(file_prefix_1, file_prefix_2);
  return file_prefix;
}

// Getter function for path to SRA FastQC analysis directory
std::pair<fs::path, fs::path> SRA::get_fastqc_dir_1() {
  std::pair<fs::path, fs::path> fastqc_dir_1(fastqc_dir_1_1, fastqc_dir_2_1);
  return fastqc_dir_1;
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

std::pair<fs::path, fs::path> SRA::get_sra_path_corr_fix() {
  std::pair<fs::path, fs::path> sra_path_corr_fix(sra_path_corr_fix_1, sra_path_corr_fix_2);
  return sra_path_corr_fix;
}

// Getter function for path to trimmed unpaired SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_trim_u() {
  std::pair<fs::path, fs::path> sra_path_trim_u(sra_path_trim_u1, sra_path_trim_u2);
  return sra_path_trim_u;
}

// Getter function for path to trimmed paired SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_trim_p() {
  std::pair<fs::path, fs::path> sra_path_trim_p(sra_path_trim_p1, sra_path_trim_p2);
  return sra_path_trim_p;
}

// Getter function for path to foreign-filtered SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_for_filt() {
  std::pair<fs::path, fs::path> sra_path_for_filt(sra_path_for_filt_1, sra_path_for_filt_2);
  return sra_path_for_filt;
}

// Getter function for fastqc directory 2
std::pair<fs::path, fs::path> SRA::get_fastqc_dir_2() {
  std::pair<fs::path, fs::path> fastqc_dir_2(fastqc_dir_1_2, fastqc_dir_2_2);
  return fastqc_dir_2;
}

// Getter function for path to overrep-filtered SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_orep_filt() {
  std::pair<fs::path, fs::path> sra_path_orep_filt(sra_path_orep_filt_1, sra_path_orep_filt_2);
  return sra_path_orep_filt;
}

// Setter function for path to SRA FastQC analysis directory
void SRA::set_fastqc_dir_1(std::pair<fs::path, fs::path> fastqcFiles1) {
  fastqc_dir_1_1 = fastqcFiles1.first;
  fastqc_dir_2_1 = fastqcFiles1.second;
}

// Setter function for path to raw SRA reads sequence
void SRA::set_sra_path_raw(std::pair<fs::path, fs::path> sraRawFiles) {
  sra_path_raw_1 = sraRawFiles.first;
  sra_path_raw_2 = sraRawFiles.second;
}

// Setter function for path to error-corrected SRA sequence
void SRA::set_sra_path_corr(std::pair<fs::path, fs::path> sraCorrFiles) {
  sra_path_corr_1 = sraCorrFiles.first;
  sra_path_corr_2 = sraCorrFiles.second;
}

// Setter function for path to unfixable error-fixed SRA sequence
void SRA::set_sra_path_corr_fix(std::pair<fs::path, fs::path> sraCorrFixFiles) {
  sra_path_corr_fix_1 = sraCorrFixFiles.first;
  sra_path_corr_fix_2 = sraCorrFixFiles.second;
}

// Setter function for path to trimmed unpaired SRA sequence
void SRA::set_sra_path_trim_u(std::pair<fs::path, fs::path> sraTrimFilesU) {
  sra_path_trim_u1 = sraTrimFilesU.first;
  sra_path_trim_u2 = sraTrimFilesU.second;
}

// Setter function for path to trimmed paired SRA sequence
void SRA::set_sra_path_trim_p(std::pair<fs::path, fs::path> sraTrimFilesP) {
  sra_path_trim_p1 = sraTrimFilesP.first;
  sra_path_trim_p2 = sraTrimFilesP.second;
}

// Setter function for foreign-filtered SRA sequence
void SRA::set_sra_path_for_filt(std::pair<fs::path, fs::path> sraForFiltFiles) {
  sra_path_for_filt_1 = sraForFiltFiles.first;
  sra_path_for_filt_2 = sraForFiltFiles.second;
}

// Setter function for fastqc directory 2
void SRA::set_fastqc_dir_2(std::pair<fs::path, fs::path> fastqcFiles2) {
  fastqc_dir_1_2 = fastqcFiles2.first;
  fastqc_dir_2_2 = fastqcFiles2.second;
}

// Setter function for path to overrep-filtered SRA sequence
void SRA::set_sra_path_orep_filt(std::pair<fs::path, fs::path> sraOrepFiltFiles) {
  sra_path_orep_filt_1 = sraOrepFiltFiles.first;
  sra_path_orep_filt_2 = sraOrepFiltFiles.second;
}

// Utility function to construct filename for SRA
std::string SRA::make_file_str() {
  std::string filename;
  filename = get_accession() + "_" + get_tax_id() + "_" +
             get_org_name().replace(get_org_name().find(" "), 1, "_");
  return filename;
}

std::string SRA::makeCheckpointName(std::string ext) {
  fs::path outDir;
  std::string cpFileName;
  if (get_accession() == "") {
    outDir = get_sra_path_corr().first.parent_path().parent_path() / "checkpoints";
    cpFileName = std::string(outDir.c_str()) + "/" + get_file_prefix().first + "." + ext + ".ok";
  }
  else {
    outDir = get_sra_path_raw().first.parent_path().parent_path() / "checkpoints";
    cpFileName = std::string(outDir.c_str()) + "/" + make_file_str() + "." + ext + ".ok";
  }
  return cpFileName;
}

void SRA::makeCheckpoint(std::string ext) {
  std::string cpFileName = makeCheckpointName(ext);
  std::ofstream cpFile;
  cpFile.open(cpFileName);
  cpFile.close();
}

bool SRA::checkpointExists(std::string ext) {
  fs::path cpFilePath(makeCheckpointName(ext).c_str());
  if (fs::exists(cpFilePath)) {
    return true;
  }
  else {
    return false;
  }
}
