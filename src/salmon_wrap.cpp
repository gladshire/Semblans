#include "salmon_wrap.h"

void combine_paired(std::vector<SRA> sras, long long int ram_b) {
  std::string inFileStr1;
  std::string outDirStr(sras[0].get_sra_path_orep_filt().first.parent_path().c_str());
  std::string outStr(outDirStr + "/" + sras[0].get_tax_id() + "_" +
                     sras[0].get_org_name().replace(sras[0].get_org_name().find(" "), 1, "_"));
  std::string outFileStr1 = outStr + "_1.fasta";
  std::string outFileStr2 = outStr + "_2.fasta";
  std::string inFileBuffer;
  inFileBuffer.reserve(ram_b);
  std::streamsize s;

  std::ofstream outFile1;
  std::ifstream inFile1;
  
  std::ofstream outFile2;
  std::ifstream inFile2;

  for (auto &sra : sras) {
    inFile1.open(sra.get_sra_path_orep_filt().first.c_str());
    outFile1.open(outFileStr1, std::ios_base::app);
    while (!inFile1.eof()) {
      inFile1.read(&inFileBuffer[0], ram_b);
      s = inFile1.gcount();
      outFile1.write(&inFileBuffer[0], s);
    }
    inFile1.close();
    if (sra.is_paired()) {
      outFile2.open(outFileStr2, std::ios_base::app);
      inFile2.open(sra.get_sra_path_orep_filt().second.c_str());
      while (!inFile2.eof()) {
        inFile2.read(&inFileBuffer[0], ram_b);
        s = inFile2.gcount();
        outFile2.write(&inFileBuffer[0], s);
      }
      inFile1.close();
      inFile2.close();
    }
    outFile1.close();
    outFile2.close();
  }
}

// Index reads from sequence data
void salmon_index(transcript trans, std::string threads,
                  bool dispOutput, std::string logFile) {
  std::string transFilePath(trans.get_trans_path_chimera().c_str());
  std::string indexFilePath(trans.get_trans_path_index().c_str());
  std::string printOut;
  if (dispOutput) {
    printOut = " |& tee -a " + logFile;
  }
  else {
    printOut = " &>> " + logFile;
  }
  int result;
  std::cout << "Generating index for transcript ..." << std::endl;

  if (fs::exists(trans.get_trans_path_index().c_str())) {
    std::cout << "Salmon index found for: " << trans.make_file_str() << std::endl;
    return;
  }
  std::string salm_cmd = PATH_SALMON + " index" + " -t " + transFilePath +
                         " -i " + indexFilePath + " -p " + threads +
                         printOut;
  result = system(salm_cmd.c_str());
  if (WIFSIGNALED(result)) {
    std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
    exit(1);
  }
}

bool runPaired(std::vector<SRA> sras) {
  int numSingle = 0;
  int numPaired = 0;
  for (auto sra : sras) {
    if (sra.is_paired()) {
      numPaired++;
    }
    else {
      numSingle++;
    }
  }
  if (numPaired >= numSingle) {
    return true;
  }
  else {
    return false;
  }
}

// Quantify reads from sequence data. Assumes all reads are either single or paired-end
void salmon_quant(transcript trans, std::vector<SRA> sras, std::string threads,
                  bool dispOutput, std::string logFile) {
  std::string transFilePath(trans.get_trans_path_chimera().c_str());
  std::string indexFilePath(trans.get_trans_path_index().c_str());
  std::string quantFilePath(trans.get_trans_path_quant().c_str());
  
  std::string sras1 = "";
  std::string sras2 = "";

  std::cout << "Generating quant for transcript ..." << std::endl;
  for (auto sraIter = sras.begin(); sraIter != sras.end(); sraIter++) {
    if (sraIter != sras.begin()) {
      sras1 += " ";
      if (sraIter->is_paired()) {
        sras2 += " ";
      }
    }
    sras1 += sraIter->get_sra_path_orep_filt().first.c_str();
    if (sraIter->is_paired()) {
      sras2 += sraIter->get_sra_path_orep_filt().second.c_str();
    }
  }
  std::string printOut;
  if (dispOutput) {
    printOut = " |& tee -a " + logFile;
  }
  else {
    printOut = " &>> " + logFile;
  }
  int result;
  if (fs::exists(trans.get_trans_path_quant())) {
    std::cout << "Quant found for: " << trans.get_org_name() << std::endl;
    return;
  }
  bool morePaired = false;
  morePaired = runPaired(sras);
  if (morePaired) {
    std::string salm_cmd = PATH_SALMON + " quant" + " -i " + indexFilePath + " --dumpEq" +
                           " --libType" + " A" + " -p " + threads +
                           " -1 " + sras1 + " -2 " + sras2 +
                           " --validateMappings" + " -o " + quantFilePath + printOut;
    result = system(salm_cmd.c_str());
  }
  else {
    std::string salm_cmd = PATH_SALMON + " quant" + " -i " + indexFilePath + " --dumpEq" +
                           " --libType" + " A" + " -p " + threads +
                           " -r " + sras1 + " --validateMappings" + " -o " + quantFilePath +
                           printOut;
    result = system(salm_cmd.c_str());
  }
  if (WIFSIGNALED(result)) {
    std::cout << "Exited with signal " << WTERMSIG(result) << std::endl;
    exit(1);
  }
}
