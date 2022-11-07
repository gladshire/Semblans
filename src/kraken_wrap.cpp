#include "kraken_wrap.h"

std::vector<std::string> get_kraken2_dbs(const INI_MAP &iniFile) {
  std::string krakDbDir = iniFile.at("Kraken2 settings").at("db_directory");
  std::vector<std::string> kraken2Dbs;
  std::string dbPath;
  for (auto db : iniFile.at("Kraken2 filter order")) {
    dbPath = krakDbDir + db.first;
    kraken2Dbs.push_back(dbPath);
  }
  return kraken2Dbs;
}


std::string get_kraken2_conf(const INI_MAP &iniFile) {
  return iniFile.at("Kraken2 settings").at("confidence_threshold");
}


void pre_summary(SRA sra, std::string db) {
  std::cout << "\nFor SRA accession: " << sra.get_accession() << std::endl;
}


void run_kraken2(std::vector<SRA> sras, std::string threads, std::string db, 
                 std::string conf_threshold, bool selfPass) {
  std::string outDir(sras[0].get_sra_path_for_filt().first.parent_path().c_str());
  std::string inFile1;
  std::string inFile2;
  std::string outFile;
  std::string repFile;
  fs::path dbPath(db.c_str());

  std::string krakCmd(PATH_KRAK2 + " --db " + db);
  std::string krakFlags("--threads " + threads + " --unclassified-out");
  for (auto sra : sras) {
    repFile = std::string(outDir + "/" + sra.make_file_str() + "." +
                          dbPath.filename().c_str() + ".report");
    if (fs::exists(fs::path(repFile.c_str()))) {
      std::cout << "Filtered version found for: " << sra.get_accession() << std::endl;
      continue;
    }
    pre_summary(sra, db);
    if (selfPass) {
      std::string tmpIn1 = std::string(sra.get_sra_path_for_filt().first.parent_path().c_str()) +
                           "/INPUT1.fq";
      std::rename(sra.get_sra_path_for_filt().first.c_str(), tmpIn1.c_str());
      inFile1 = tmpIn1;
    }
    else {
      inFile1 = sra.get_sra_path_trim_p().first.c_str();
    }
    outFile = sra.get_sra_path_for_filt().first.c_str();
    if (sra.is_paired()) {
      if (selfPass) {
        std::string tmpIn2 = std::string(sra.get_sra_path_for_filt().second.parent_path().c_str()) +
                             "/INPUT2.fq";
        std::rename(sra.get_sra_path_for_filt().second.c_str(), tmpIn2.c_str());
        inFile2 = tmpIn2;
      }
      else {
        inFile2 = sra.get_sra_path_trim_p().second.c_str();
      }
      outFile = std::string(outFile).replace(outFile.length() - 10, 2, "#");
      system((krakCmd + " " + krakFlags + " " + outFile + " --paired " + "--output - " +
              inFile1 + " " + inFile2 + " --confidence " + conf_threshold + " --report " +
              outDir + "/" + sra.make_file_str() + "." + dbPath.filename().c_str() +
              ".report").c_str());
    }
    else {
      system((krakCmd + " " + krakFlags + " " + outFile + "--output - " + inFile1 + " " +
              " --confidence " + conf_threshold + " --report " + outDir + "/" +
              sra.make_file_str() + "." + dbPath.filename().c_str() + ".report").c_str());
    }
  }
  std::string tmpIn1 = std::string(sras[0].get_sra_path_for_filt().first.parent_path().c_str()) + "/INPUT1.fq";
  std::string tmpIn2 = std::string(sras[0].get_sra_path_for_filt().first.parent_path().c_str()) + "/INPUT2.fq";
  std::remove(tmpIn1.c_str());
  std::remove(tmpIn2.c_str());
}


void run_kraken2_dbs(std::vector<SRA> sras, std::string threads, std::vector<std::string> dbs,
                     std::string conf_threshold) {
  std::cout << "\nFiltering foreign reads for:\n" << std::endl;
  summarize_all_sras(sras);
  int i = 0;
  for (auto db : dbs) {
    std::cout << "\nNow filtering: " << fs::path(db.c_str()).filename().c_str() << std::endl;
    if (i == 0) {
      run_kraken2(sras, threads, db, conf_threshold, false);
    }
    else {
      run_kraken2(sras, threads, db, conf_threshold, true);
    }
    i++;
  }
}
