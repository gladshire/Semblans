#include "sra.h"


extern std::vector<std::string> stepDirs;

// Default constructor for SRA object
SRA::SRA() {
  sra_accession = "";
  org_name = "";
  tax_id = "";
  spots = 0;
  spots_m = -1;
  bp = 0;
  paired = false;
  compressedFiles = true;
  file_prefix_1 = "";
  file_prefix_2 = "";
}

// SRA accession number constructor for SRA object
// Takes SRA accession number string as input
// Fills object members with correct values using NCBI eutils API
SRA::SRA(std::string sra_accession, INI_MAP cfgIni, bool dispOutput,
         bool compressedFiles, std::string logFile, int num) {
  std::string outDir(fs::canonical(fs::path(cfgIni["General"]["output_directory"].c_str())).c_str());
  std::string projName(cfgIni["General"]["project_name"]);
  std::string apiKey(cfgIni["General"]["ncbi_api_key"]);
  this->compressedFiles = compressedFiles;
  std::string compressExt;
  std::chrono::milliseconds queryLim(500);
  if (apiKey != "") {
    std::chrono::milliseconds queryLim(250);
  }

  // Download temp XML file for SRA accession, containing information for object members
  std::string curlCmdStr = "curl -s -o .tmp" + std::to_string(num) + ".xml \"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?rettype=runinfo&db=sra&id=\"" +
                           sra_accession + "&api_key=" + apiKey;
  int result;
  int numRetries = 0;
  std::string platform;

  logOutput("\nObtaining information for accession(s)", logFile);
  while (true) {
    // After 5 unsuccessful retries, break and exit
    if (numRetries == 5) {
      logOutput("\n  ERROR: Could not retrieve accession \"" + sra_accession + "\"\n", logFile);
      this->sra_accession = "FAILURE";
      break;
    }
    system(curlCmdStr.c_str());
    std::this_thread::sleep_for(queryLim);

    // Parse XML file for object information
    try {
      rapidxml::file<> sra_xml((".tmp" + std::to_string(num) + ".xml").c_str());
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
      if (!strcmp(parse_node->value(), "PAIRED") && spots_m == 0) {
        // NO single spots
        paired = false;
      }
      else if (!strcmp(parse_node->value(), "PAIRED") && spots_m != spots_m) {
        // There are single spots
        paired = true;
      }
      else {
        // NO single spots
        paired = true;
      }

      // Ensure sequencing platform is ILLUMINA
      while (strcmp(parse_node->name(), "Platform")) {
        parse_node = parse_node->next_sibling();
      }
      platform = std::string(parse_node->value());
      if (platform != "ILLUMINA") {
        std::cout << "ERROR: Accession \'" << sra_accession << "\' was not run on Illumina" << std::endl;
        std::cout << "  All SRA accessions must be Illumina runs, and cannot be PacBio, 10X, etc" << std::endl;
        exit(1);
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
      numRetries++;
      queryLim *= 2;
    }
    catch (const rapidxml::parse_error & e) {
      numRetries++;
      queryLim *= 2;
    }
  }
  if (this->sra_accession != "FAILURE") {
    if (dispOutput) {
      logOutput("\n  Successfully retrieved information for: " + sra_accession, logFile);
    }

    std::string fileBase = make_file_str();
    if (paired) {
      file_prefix_1 = fileBase + "_1";
      file_prefix_2 = fileBase + "_2";
    }
    else {
      file_prefix_1 = fileBase;
      file_prefix_2 = "";
    }

    std::string projPath((fs::path(outDir.c_str()) / fs::path(projName.c_str())).c_str());
    projPath += "/";

    extern std::vector<std::string> stepDirs;

    if (compressedFiles) {
      compressExt = ".gz";
    }
    else {
      compressExt = "";
    }

    sra_path_raw_1  = (projPath + stepDirs[0] + fileBase + ".fastq" + compressExt).c_str();
    fastqc_dir_1_1  = (projPath + stepDirs[1] + fileBase + "/" + fileBase).c_str();
    sra_path_corr_1 = (projPath + stepDirs[2] + fileBase + ".cor.fq" + compressExt).c_str();
    sra_path_corr_fix_1 = (projPath + stepDirs[2] + fileBase + ".cor.fix.fq" + compressExt).c_str();
    sra_path_trim_u1 = (projPath + stepDirs[3] + fileBase + ".trim.fq" + compressExt).c_str();
    sra_path_trim_p1 = (projPath + stepDirs[3] + fileBase + ".trim.fq" + compressExt).c_str();
    sra_path_for_filt_1 = (projPath + stepDirs[4] + fileBase + ".filt.fq" + compressExt).c_str();
    fastqc_dir_1_2 = (projPath + stepDirs[5] + fileBase + "/" + fileBase).c_str();
    sra_path_orep_filt_1 = (projPath + stepDirs[6] + fileBase + ".orep.filt.fq" + compressExt).c_str();

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

      int extIndMod;
      if (compressedFiles) {
	extIndMod = 3;
      }
      else {
	extIndMod = 0;
      }
      sra_path_raw_1 = (projPath + stepDirs[0] + fileBase + "_1.fastq" + compressExt).c_str();
      fastqc_dir_1_1 = (fastqc_dir_1_1_str.insert(fastqc_dir_1_1_str.length(), "_1")).c_str();
      sra_path_corr_1 = (sra_path_corr_1_str.insert(sra_path_corr_1_str.length() - 7 - extIndMod, "_1")).c_str();
      sra_path_corr_fix_1 = (sra_path_corr_fix_1_str.insert(sra_path_corr_fix_1_str.length() - 11 - extIndMod, "_1")).c_str();
      sra_path_trim_u1 = (sra_path_trim_u1_str.insert(sra_path_trim_u1_str.length() - 8 - extIndMod, "_1.unpaired")).c_str();
      sra_path_trim_p1 = (sra_path_trim_p1_str.insert(sra_path_trim_p1_str.length() - 8 - extIndMod, "_1.paired")).c_str();
      sra_path_for_filt_1 = (sra_path_for_filt_1_str.insert(sra_path_for_filt_1_str.length() - 8 - extIndMod, "_1")).c_str();
      fastqc_dir_1_2 = (fastqc_dir_1_2_str.insert(fastqc_dir_1_2_str.length(), "_1")).c_str();
      sra_path_orep_filt_1 = (sra_path_orep_filt_1_str.insert(sra_path_orep_filt_1_str.length() - 13 - extIndMod, "_1")).c_str();

      sra_path_raw_2 = (projPath + stepDirs[0] + fileBase + "_2.fastq" + compressExt).c_str();
      fastqc_dir_2_1 = (projPath + stepDirs[1] + fileBase + "/" + fileBase + "_2").c_str();
      sra_path_corr_2 = (projPath + stepDirs[2] + fileBase + "_2.cor.fq" + compressExt).c_str();
      sra_path_corr_fix_2 = (projPath + stepDirs[2] + fileBase + "_2.cor.fix.fq" + compressExt).c_str();
      sra_path_trim_u2 = (projPath + stepDirs[3] + fileBase + "_2.unpaired.trim.fq" + compressExt).c_str();
      sra_path_trim_p2 = (projPath + stepDirs[3] + fileBase + "_2.paired.trim.fq" + compressExt).c_str();
      sra_path_for_filt_2 = (projPath + stepDirs[4] + fileBase + "_2.filt.fq" + compressExt).c_str();
      fastqc_dir_2_2 = (projPath + stepDirs[5] + fileBase + "/" + fileBase + "_2").c_str();
      sra_path_orep_filt_2 = (projPath + stepDirs[6] + fileBase + "_2.orep.filt.fq" + compressExt).c_str();
    }
  }
}


// Constructor for if local paired-end data used
SRA::SRA(std::string fileName1, std::string fileName2, INI_MAP cfgIni, bool compressedFiles,
         std::string logFile) {


  std::string outDir(fs::canonical(fs::path(cfgIni["General"]["output_directory"].c_str())).c_str());
  std::string projName(cfgIni["General"]["project_name"]);
  std::string projPath((fs::path(outDir.c_str()) / fs::path(projName.c_str())).c_str());
  projPath += "/";

  std::string localDataDir(cfgIni["General"]["local_data_directory"]);
  std::string compressExt;
  std::string fileBase1(fs::path(fileName1.c_str()).stem().c_str());
  std::string fileBase2(fs::path(fileName2.c_str()).stem().c_str());

  file_prefix_1 = fileBase1;
  file_prefix_2 = fileBase2;

  sra_accession = "";
  org_name = "";
  tax_id = "";
  spots = 0;
  spots_m = -1;
  bp = 0;
  paired = (fileName2 == "") ? false : true;
  this->compressedFiles = compressedFiles;

  if (compressedFiles) {
    compressExt = ".gz";
  }
  else {
    compressExt = "";
  }

  std::ifstream sraFile1;
  std::ifstream sraFile2;
  std::streamsize s;

  uintmax_t numReads1 = 0;
  uintmax_t numReads2 = 0;

  std::string buffer;
  buffer.reserve(1000000000);

  sraFile1.open(fileName1);
  if (paired) {
    sraFile2.open(fileName2);
  }

  char * nlPos;
  char * inFileL;

  while (!sraFile1.eof() && !sraFile1.good()) {
    sraFile1.read(&buffer[0], 1000000000);
    s = sraFile1.gcount();

    nlPos = &buffer[0];
    inFileL = &buffer[0] + s;
    while (nlPos != inFileL) {
      nlPos = std::find(nlPos + 1, inFileL, '\n');
      if (nlPos != inFileL) {
        numReads1++;
      }
    }
  }
  numReads1 /= 4;

  if (paired) {
    while (!sraFile2.eof() && !sraFile2.good()) {
      sraFile2.read(&buffer[0], 1000000000);
      s = sraFile2.gcount();

      nlPos = &buffer[0];
      inFileL = &buffer[0] + s;
      while (nlPos != inFileL) {
        nlPos = std::find(nlPos + 1, inFileL, '\n');
        if (nlPos != inFileL) {
          numReads2++;
        }
      }
    }
    numReads2 /= 4;
  }

  if (paired && (numReads1 != numReads2)) {
    logOutput("ERROR: Forward and reverse files do not have the same number of reads.", logFile);
    exit(1);
  }

  spots = numReads1;

  sra_path_raw_1 = fileName1.c_str();
  fastqc_dir_1_1  = (projPath + stepDirs[1] + fileBase1 + "/" + fileBase1).c_str();
  sra_path_corr_1 = (projPath + stepDirs[2] + fileBase1 + ".cor.fq" + compressExt).c_str();
  sra_path_corr_fix_1 = (projPath + stepDirs[2] + fileBase1 + ".cor.fix.fq" + compressExt).c_str();
  sra_path_trim_u1 = (projPath + stepDirs[3] + fileBase1 + ".trim.fq" + compressExt).c_str();
  sra_path_trim_p1 = (projPath + stepDirs[3] + fileBase1 + ".trim.fq" + compressExt).c_str();
  sra_path_for_filt_1 = (projPath + stepDirs[4] + fileBase1 + ".filt.fq" + compressExt).c_str();
  fastqc_dir_1_2 = (projPath + stepDirs[5] + fileBase1 + "/" + fileBase1).c_str();
  sra_path_orep_filt_1 = (projPath + stepDirs[6] + fileBase1 + ".orep.filt.fq" + compressExt).c_str();

  if (paired) {
    sra_path_raw_2 = fileName2.c_str();
    fastqc_dir_2_1 = (projPath + stepDirs[1] + fileBase2 + "/" + fileBase2).c_str();
    sra_path_corr_2 = (projPath + stepDirs[2] + fileBase2 + ".cor.fq" + compressExt).c_str();
    sra_path_corr_fix_2 = (projPath + stepDirs[2] + fileBase2 + ".cor.fix.fq" + compressExt).c_str();
    sra_path_trim_u1 = (projPath + stepDirs[3] + fileBase1 + ".unpaired.trim.fq" + compressExt).c_str();
    sra_path_trim_p1 = (projPath + stepDirs[3] + fileBase1 + ".paired.trim.fq" + compressExt).c_str();
    sra_path_trim_u2 = (projPath + stepDirs[3] + fileBase2 + ".unpaired.trim.fq" + compressExt).c_str();
    sra_path_trim_p2 = (projPath + stepDirs[3] + fileBase2 + ".paired.trim.fq" + compressExt).c_str();
    sra_path_for_filt_2 = (projPath + stepDirs[4] + fileBase2 + ".filt.fq" + compressExt).c_str();
    fastqc_dir_2_2 = (projPath + stepDirs[5] + fileBase2 + "/" + fileBase2).c_str();
    sra_path_orep_filt_2 = (projPath + stepDirs[6] + fileBase2 + ".orep.filt.fq" + compressExt).c_str();
  }
}

SRA::SRA(std::string fileName1, std::string fileName2, std::string outDir, bool compressedFiles,
         bool dummy) {

  sra_accession = "";
  org_name = "";
  tax_id = "";
  spots = 0;
  spots_m = -1;
  bp = 0;
  paired = (fileName2 == "null") ? false : true;
  this->compressedFiles = compressedFiles;

  //outDir = std::string(fs::canonical(outDir.c_str()).parent_path().c_str()) + "/";
  fileName1 = std::string(fs::canonical(fileName1.c_str()).c_str());
  if (paired) {
    fileName2 = std::string(fs::canonical(fileName2.c_str()).c_str());
  }

  //std::string localDataDir(cfgIni["General"]["local_data_directory"]);
  std::string compressExt;
  std::string fileBase1(fs::path(fileName1.c_str()).stem().c_str());
  std::string fileBase2(fs::path(fileName2.c_str()).stem().c_str());

  file_prefix_1 = fileBase1;
  file_prefix_2 = (paired) ? fileBase2 : "";


  if (compressedFiles) {
    compressExt = ".gz";
  }
  else {
    compressExt = "";
  }

  std::ifstream sraFile1;
  std::ifstream sraFile2;
  std::streamsize s;

  uintmax_t numReads1 = 0;
  uintmax_t numReads2 = 0;

  std::string buffer;
  buffer.reserve(1000000000);

  sraFile1.open(fileName1);
  if (paired) {
    sraFile2.open(fileName2);
  }

  char * nlPos;
  char * inFileL;

  while (!sraFile1.eof() && !sraFile1.good()) {
    sraFile1.read(&buffer[0], 1000000000);
    s = sraFile1.gcount();

    nlPos = &buffer[0];
    inFileL = &buffer[0] + s;
    while (nlPos != inFileL) {
      nlPos = std::find(nlPos + 1, inFileL, '\n');
      if (nlPos != inFileL) {
        numReads1++;
      }
    }
  }
  numReads1 /= 4;

  if (paired) {
    while (!sraFile2.eof() && !sraFile2.good()) {
      sraFile2.read(&buffer[0], 1000000000);
      s = sraFile2.gcount();

      nlPos = &buffer[0];
      inFileL = &buffer[0] + s;
      while (nlPos != inFileL) {
        nlPos = std::find(nlPos + 1, inFileL, '\n');
        if (nlPos != inFileL) {
          numReads2++;
        }
      }
    }
    numReads2 /= 4;
  }

  spots = numReads1;

  sra_path_raw_1 = fileName1.c_str();
  if (!dummy) {
    fastqc_dir_1_1  = (outDir + stepDirs[1] + fileBase1 + "/" + fileBase1).c_str();
    sra_path_corr_1 = (outDir + stepDirs[2] + fileBase1 + ".cor.fq" + compressExt).c_str();
    sra_path_corr_fix_1 = (outDir + stepDirs[2] + fileBase1 + ".cor.fix.fq" + compressExt).c_str();
    sra_path_trim_u1 = (outDir + stepDirs[3] + fileBase1 + ".trim.fq" + compressExt).c_str();
    sra_path_trim_p1 = (outDir + stepDirs[3] + fileBase1 + ".trim.fq" + compressExt).c_str();
    sra_path_for_filt_1 = (outDir + stepDirs[4] + fileBase1 + ".filt.fq" + compressExt).c_str();
    fastqc_dir_1_2 = (outDir + stepDirs[5] + fileBase1 + "/" + fileBase1).c_str();
    sra_path_orep_filt_1 = (outDir + stepDirs[6] + fileBase1 + ".orep.filt.fq" + compressExt).c_str();
  }
  if (paired) {
    sra_path_raw_2 = fileName2.c_str();
    if (!dummy) {
      fastqc_dir_2_1 = (outDir + stepDirs[1] + fileBase2 + "/" + fileBase2).c_str();
      sra_path_corr_2 = (outDir + stepDirs[2] + fileBase2 + ".cor.fq" + compressExt).c_str();
      sra_path_corr_fix_2 = (outDir + stepDirs[2] + fileBase2 + ".cor.fix.fq" + compressExt).c_str();
      sra_path_trim_u1 = (outDir + stepDirs[3] + fileBase1 + ".unpaired.trim.fq" + compressExt).c_str();
      sra_path_trim_p1 = (outDir + stepDirs[3] + fileBase1 + ".paired.trim.fq" + compressExt).c_str();
      sra_path_trim_u2 = (outDir + stepDirs[3] + fileBase2 + ".unpaired.trim.fq" + compressExt).c_str();
      sra_path_trim_p2 = (outDir + stepDirs[3] + fileBase2 + ".paired.trim.fq" + compressExt).c_str();
      sra_path_for_filt_2 = (outDir + stepDirs[4] + fileBase2 + ".filt.fq" + compressExt).c_str();
      fastqc_dir_2_2 = (outDir + stepDirs[5] + fileBase2 + "/" + fileBase2).c_str();
      sra_path_orep_filt_2 = (outDir + stepDirs[6] + fileBase2 + ".orep.filt.fq" + compressExt).c_str();
    }
  }
}

// Copy constructor for SRA object
SRA::SRA(const SRA & sra) {
  sra_accession = sra.get_accession();
  org_name = sra.get_org_name();
  tax_id = sra.get_tax_id();
  spots = sra.get_num_reads();
  spots_m = sra.get_bp();
  paired = sra.is_paired();
  compressedFiles = false;
  file_prefix_1 = sra.get_file_prefix().first;
  file_prefix_2 = sra.get_file_prefix().second;
  fastqc_dir_1_1 = sra.get_fastqc_dir_1().first;
  fastqc_dir_2_1 = sra.get_fastqc_dir_1().second;
  sra_path_raw_1 = sra.get_sra_path_raw().first;
  sra_path_raw_2 = sra.get_sra_path_raw().second;
  sra_path_corr_1 = sra.get_sra_path_corr().first;
  sra_path_corr_2 = sra.get_sra_path_corr().second;
  sra_path_corr_fix_1 = sra.get_sra_path_corr_fix().first;
  sra_path_corr_fix_2 = sra.get_sra_path_corr_fix().second;
  sra_path_trim_u1 = sra.get_sra_path_trim_u().first;
  sra_path_trim_u2 = sra.get_sra_path_trim_u().second;
  sra_path_trim_p1 = sra.get_sra_path_trim_p().first;
  sra_path_trim_p2 = sra.get_sra_path_trim_p().second;
  sra_path_for_filt_1 = sra.get_sra_path_for_filt().first;
  sra_path_for_filt_2 = sra.get_sra_path_for_filt().second;
  fastqc_dir_1_2 = sra.get_fastqc_dir_2().first;
  fastqc_dir_2_2 = sra.get_fastqc_dir_2().second;
  sra_path_orep_filt_1 = sra.get_sra_path_orep_filt().first;
  sra_path_orep_filt_2 = sra.get_sra_path_orep_filt().second;
}

// Getter function for object SRA accession number
std::string SRA::get_accession() const {
  return sra_accession;
}

// Getter function for object organism name
std::string SRA::get_org_name() const {
  return org_name;
}

// Getter function for object taxonomic ID
std::string SRA::get_tax_id() const {
  return tax_id;
}

// Getter function for object number of read spots
long unsigned int SRA::get_num_reads() const {
  return spots;
}

void SRA::set_num_reads(long unsigned int numReads) {
  spots = numReads;
}

// Getter function for object number of read spots that are paired
int SRA::get_spots_m() const {
  return spots_m;
}

// Getter function for object number of base pairs in SRA run
uintmax_t SRA::get_bp() const {
  return bp;
}

// Getter function for object single-end or paired-end layout
bool SRA::is_paired() const {
  return paired;
}

std::pair<std::string, std::string> SRA::get_file_prefix() const {
  std::pair<std::string, std::string> file_prefix(file_prefix_1, file_prefix_2);
  return file_prefix;
}

// Getter function for path to SRA FastQC analysis directory
std::pair<fs::path, fs::path> SRA::get_fastqc_dir_1() const {
  std::pair<fs::path, fs::path> fastqc_dir_1(fastqc_dir_1_1, fastqc_dir_2_1);
  return fastqc_dir_1;
}

// Getter function for path to raw SRA reads sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_raw() const {
  std::pair<fs::path, fs::path> sra_path_raw(sra_path_raw_1, sra_path_raw_2);
  return sra_path_raw;
}

// Getter function for path to error-corrected SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_corr() const {
  std::pair<fs::path, fs::path> sra_path_corr(sra_path_corr_1, sra_path_corr_2);
  return sra_path_corr;
}

std::pair<fs::path, fs::path> SRA::get_sra_path_corr_fix() const {
  std::pair<fs::path, fs::path> sra_path_corr_fix(sra_path_corr_fix_1, sra_path_corr_fix_2);
  return sra_path_corr_fix;
}

// Getter function for path to trimmed unpaired SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_trim_u() const {
  std::pair<fs::path, fs::path> sra_path_trim_u(sra_path_trim_u1, sra_path_trim_u2);
  return sra_path_trim_u;
}

// Getter function for path to trimmed paired SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_trim_p() const {
  std::pair<fs::path, fs::path> sra_path_trim_p(sra_path_trim_p1, sra_path_trim_p2);
  return sra_path_trim_p;
}

// Getter function for path to foreign-filtered SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_for_filt() const {
  std::pair<fs::path, fs::path> sra_path_for_filt(sra_path_for_filt_1, sra_path_for_filt_2);
  return sra_path_for_filt;
}

// Getter function for fastqc directory 2
std::pair<fs::path, fs::path> SRA::get_fastqc_dir_2() const {
  std::pair<fs::path, fs::path> fastqc_dir_2(fastqc_dir_1_2, fastqc_dir_2_2);
  return fastqc_dir_2;
}

// Getter function for path to overrep-filtered SRA sequence
std::pair<fs::path, fs::path> SRA::get_sra_path_orep_filt() const {
  std::pair<fs::path, fs::path> sra_path_orep_filt(sra_path_orep_filt_1, sra_path_orep_filt_2);
  return sra_path_orep_filt;
}

// Setter function for number of reads in SRA run
void SRA::set_num_reads(long int numReads) {
  this->spots = numReads;
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
  std::string filenameCleaned;
  std::string orgName = get_org_name();
  while (orgName.find(" ") != std::string::npos) {
    orgName = orgName.replace(orgName.find(" "), 1, "_");
  }
  filename = get_accession() + "_" + get_tax_id() + "_" + orgName;
  for (int i = 0; i < filename.size(); i++) {
    if (isalnum(filename[i]) || filename[i] == '_') {
      filenameCleaned.push_back(filename[i]);
    }
  }
  return filenameCleaned;
}

std::string SRA::makeCheckpointName(std::string ext) {
  fs::path outDir;
  std::string cpFileName;
  outDir = get_fastqc_dir_2().first.parent_path().parent_path().parent_path().parent_path() / ".checkpoints";
  if (get_accession() == "") {
    cpFileName = std::string(outDir.c_str()) + "/" + get_file_prefix().first + "." + ext + ".ok";
  }
  else {
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
