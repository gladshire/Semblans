#include "kraken_wrap.h"

void run_kraken2(std::vector<SRA> sras, std::string threads, std::string db, bool selfPass) {
  std::string outDir(sras[0].get_sra_path_filt().first.parent_path().c_str());
  std::string inFile1;
  std::string inFile2;
  std::string outFile;
  fs::path dbPath(db.c_str());

  std::string krakCmd(PATH_KRAK2 + " --db " + db);
  std::string krakFlags("--threads " + threads + " --unclassified-out");
  for (auto sra : sras) {
    if (selfPass) {
      std::string tmpIn1 = std::string(sra.get_sra_path_filt().first.parent_path().c_str()) +
                           "/INPUT1.fq";
      std::rename(sra.get_sra_path_filt().first.c_str(), tmpIn1.c_str());
      inFile1 = tmpIn1;
    }
    else {
      inFile1 = sra.get_sra_path_trim_p().first.c_str();
    }
    outFile = sra.get_sra_path_filt().first.c_str();
    if (sra.is_paired()) {
      if (selfPass) {
        std::string tmpIn2 = std::string(sra.get_sra_path_filt().second.parent_path().c_str()) +
                             "/INPUT2.fq";
        std::rename(sra.get_sra_path_filt().second.c_str(), tmpIn2.c_str());
        inFile2 = tmpIn2;
      }
      else {
        inFile2 = sra.get_sra_path_trim_p().second.c_str();
      }
      outFile = std::string(outFile).replace(outFile.length() - 10, 2, "#");
      if (fs::exists(sra.get_sra_path_filt().first) &&
          fs::exists(sra.get_sra_path_filt().second)) {
        std::cout << "Filtered version found for: " << sra.get_accession() << std::endl;
        continue;
      }
      system((krakCmd + " " + krakFlags + " " + outFile + " --paired" + "--output - " +
              inFile1 + " " + inFile2 + " --report " + outDir + "/" + sra.make_file_str() +
              "." + dbPath.filename().c_str() + ".report").c_str());
    }
    else {
      if (fs::exists(sra.get_sra_path_filt().first)) {
        std::cout << "Filtered version found for: " << sra.get_accession() << std::endl;
        continue;
      }
      system((krakCmd + " " + krakFlags + " " + outFile + "--output - " + inFile1 + " " +
              " --report " + outDir + "/" + sra.make_file_str() + "." +
              dbPath.filename().c_str() + ".report").c_str());
    }
  }
}


void run_kraken2_dbs(std::vector<SRA> sras, std::string threads, std::vector<std::string> dbs) {
  int i = 0;
  for (auto db : dbs) {
    if (i == 0) {
      run_kraken2(sras, threads, db, false);
    }
    else {
      run_kraken2(sras, threads, db, true);
    }
  }
}
