#!/bin/bash

# rcorrVersion   =
# trimmVersion   =
# fastqcVersion  = 
# kraken2Version =
# trinityVersion = 
# dmndVersion    = 
# blastVersion   =
# corsetVersion  =
# salmonVersion  = 
# trnsDecVersion =

echo "Initiating install of Semblans"

# Prepare Semblans file structure
mkdir include
mkdir lib
mkdir external
mkdir data

#========================================================================
#/////////////   INSTALLATION OF THIRD-PARTY LIBRARIES   \\\\\\\\\\\\\\\\
#========================================================================

echo "Now installing required libraries"

# Install boost libraries
echo "  Installing Boost libraries ..."
wget -q https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
tar -xf boost_1_81_0.tar.gz
cd boost_1_81_0
./bootstrap.sh --prefix=../
./b2 install
mv LICENSE_1_0.txt ../include/boost/
cd ..
rm -rf boost_1_81_0*

# Install rapidXML
echo "  Installing rapidXML library ..."
wget -q https://sourceforge.net/projects/rapidxml/files/rapidxml/rapidxml%201.13/rapidxml-1.13.zip
unzip rapidxml-1.13.zip -d ./include/
mv ./include/rapidxml-1.13 ./include/rapidxml
rm rapidxml-1.13.zip

# Install libconfini
echo "  Installing libconfini library ..."
wget -q https://github.com/madmurphy/libconfini/releases/download/1.16.4/libconfini-1.16.4-x86_64-bin.tar.xz
tar -xf libconfini-1.16.4-x86_64-bin.tar.xz
mkdir ./include/libconfini
mv ./usr/include/* ./include/libconfini/
mv ./usr/lib/* ./lib/
mv ./usr/share/doc/libconfini/AUTHORS ./include/libconfini/
mv ./usr/share/doc/libconfini/COPYING ./include/libconfini/
rm libconfini-1.16.4-x86_64-bin.tar.xz
rm -rf ./usr/

# Install libcurl
echo " Installing libcurl library ..."
wget -q https://curl.se/download/curl-8.1.2.tar.gz
tar -xf curl-8.1.2.tar.gz
mv ./curl-8.1.2/include/curl ./include/
mv ./curl-8.1.2/COPYING ./include/curl
rm -rf curl-8.1.2
rm curl-8.1.2.tar.gz

#========================================================================
#/////////////   INSTALLATION OF THIRD-PARTY PACKAGES    \\\\\\\\\\\\\\\\
#========================================================================

echo "Now installing required packages ..."

# Install HMMER
echo "Installing HMMER ..."
wget -q http://eddylab.org/software/hmmer/hmmer.tar.gz
tar -xf hmmer.tar.gz -C ./external/
mv ./external/hmmer* ./external/hmmer
cd ./external/hmmer
./configure --prefix "$PWD"
make
make install
cd ../../
rm hmmer.tar.gz

# Install PANTHER HMM Scoring components
echo "Downloading PANTHER components ..."
wget -q http://data.pantherdb.org/ftp/panther_library/current_release/PANTHER17.0_hmmscoring.tgz
tar -xzf PANTHER17.0_hmmscoring.tgz
mv target panther_db/
mv panther_db ./data/
rm PANTHER17.0_hmmscoring.tgz
wget -r -q --no-parent http://data.pantherdb.org/ftp/hmm_scoring/current_release/pantherScore2.2/
mv data.pantherdb.org/ftp/hmm_scoring/current_release/pantherScore2.2 ./external/pantherScore
chmod +x ./external/pantherScore/pantherScore2.2.pl
rm -rf data.pantherdb.org


# Install NCBI sra-tools
echo "Installing SRA-Tools ..."
wget -q --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xf sratoolkit.tar.gz -C ./external/
mv ./external/sratoolkit* ./external/sratoolkit
cd ./external/sratoolkit
wget -q https://raw.githubusercontent.com/ncbi/sra-tools/master/LICENSE
cd ../../
rm sratoolkit.tar.gz

# Install FastQC
echo "Installing FastQC ..."
wget -q --output-document fastqc.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc.zip -d ./external/
chmod 755 ./external/FastQC/fastqc
rm fastqc.zip

# Install Rcorrector
echo "Installing Rcorrector ..."
wget -q --output-document rcorrector.tar.gz https://github.com/mourisl/Rcorrector/archive/refs/tags/v1.0.5.tar.gz
tar -xf rcorrector.tar.gz -C ./external/
mv ./external/Rcorrector-1.0.5 ./external/Rcorrector
cd ./external/Rcorrector
make
cd ../../
rm rcorrector.tar.gz

# Install Trimmomatic
echo "Installing Trimmomatic ..."
wget -q --output-document trimmomatic.zip https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip
unzip trimmomatic.zip -d ./external/
mv ./external/Trimmomatic-0.39 ./external/Trimmomatic
cd ./external/Trimmomatic/adapters/
echo -e "$(cat NexteraPE-PE.fa)\n$(cat TruSeq2-PE.fa TruSeq2-SE.fa TruSeq3-PE-2.fa)\n$(cat TruSeq3-PE.fa)\n$(cat TruSeq3-SE.fa)" > TruSeq_all.fa
cd ../../../
rm trimmomatic.zip

# Install Kraken2
echo "Installing Kraken2 ..."
wget -q --output-document kraken2.tar.gz https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.tar.gz
tar -xf kraken2.tar.gz -C ./external/
mv ./external/kraken2-2.1.2 ./external/kraken2
cd ./external/kraken2
./install_kraken2.sh .
cd ../../
rm kraken2.tar.gz

# Install Trinity
echo "Installing Trinity ..."
wget -q --output-document trinity.tar.gz https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz
tar -xf trinity.tar.gz -C ./external/
mv ./external/trinityrnaseq-v2.15.1 ./external/trinityrnaseq
cd ./external/trinityrnaseq/trinity-plugins/bamsifter
sed -i '1s/^/#include <string> \n/' sift_bam_max_cov.cpp
cd ../../
make
cd ../../
rm trinity.tar.gz

# Install NCBI BLAST
echo "Installing NCBI-BLAST+ ..."
wget -q --output-document ncbi-blast.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.1+-x64-linux.tar.gz
tar -xf ncbi-blast.tar.gz -C ./external/
mv ./external/ncbi-blast* ./external/ncbi-blast+
rm ncbi-blast.tar.gz

# Install Diamond
echo "Installing Diamond ..."
mkdir ./external/diamond-2.1.7
wget -q --output-document diamond.tar.gz https://github.com/bbuchfink/diamond/releases/download/v2.1.7/diamond-linux64.tar.gz
tar -xf diamond.tar.gz -C ./external/diamond-2.1.7/
mv ./external/diamond-2.1.7 ./external/diamond
cd ./external/diamond
wget -q https://raw.githubusercontent.com/bbuchfink/diamond/master/LICENSE
cd ../../
rm diamond.tar.gz


# Install Corset
echo "Installing Corset ..."
wget -q --output-document corset.tar.gz https://github.com/Oshlack/Corset/releases/download/version-1.09/corset-1.09-linux64.tar.gz
tar -xf corset.tar.gz -C ./external/
mv ./external/corset-1.09-linux64 ./external/corset
cd ./external/corset/
chmod -x LICENSE
chmod -x COPYING
cd ../../
rm corset.tar.gz

# Install STAR
echo "Installing STAR ..."
wget -q --output-document star.tar.gz https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11a.tar.gz
tar -xf star.tar.gz -C ./external/
mv ./external/STAR* ./external/STAR
cd ./external/STAR/source/
make STAR
cd ../../../
rm star.tar.gz

# Install Salmon
echo "Installing Salmon ..."
wget -q --output-document salmon.tar.gz https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
tar -xf salmon.tar.gz -C ./external/
mv ./external/salmon* ./external/salmon
cd ./external/salmon
wget -q https://raw.githubusercontent.com/COMBINE-lab/salmon/master/LICENSE
cd ../../
rm salmon.tar.gz

# Install TransDecoder
echo "Installing TransDecoder ..."
wget -q --output-document transdecoder.tar.gz https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.0.tar.gz
tar -xf transdecoder.tar.gz -C ./external/
mv ./external/TransDecoder-TransDecoder-v5.7.0 ./external/TransDecoder
rm transdecoder.tar.gz

#========================================================================
#//////////   TEST SUCCESSFUL INSTALLATION OF DEPENDENCIES    \\\\\\\\\\\
#========================================================================

packagesNotInstalled=()

# Check HMMER installation
if [ ! -e "./external/hmmer/bin/hmmscan" ]; then
	packagesNotInstalled+=("HMMER")
fi

# Check pantherScore installation
if [ ! -e "./external/pantherScore/pantherScore2.2.pl" ]; then
	packagesNotInstalled+=("pantherScore")
fi

# Check SRA Toolkit installation
if ( [ ! -e "./external/sratoolkit/bin/prefetch" ] &&
     [ -e "./external/sratoolkit/bin/fasterq-dump" ] ); then
	packagesNotInstalled+=("SRA-Tools")
fi

# Check FastQC installation
if [ ! -e "./external/FastQC/fastqc" ]; then
	packagesNotInstalled+=("FastQC")
fi

# Check Rcorrector installation
if [ ! -e "./external/Rcorrector/rcorrector" ]; then
	packagesNotInstalled+=("Rcorrector")
fi

# Check Trimmomatic installation
if [ ! -e "./external/Trimmomatic/trimmomatic-0.39.jar" ]; then
	packagesNotInstalled+=("Trimmomatic")
fi

# Check Kraken2 installation
if [ ! -e "./external/kraken2/kraken2" ]; then
	packagesNotInstalled+=("Kraken2")
fi

# Check Trinity installation
if [ ! -e "./external/trinityrnaseq/Trinity" ]; then
	packagesNotInstalled+=("Trinity")
fi

# Check BLAST installation
if ( [ ! -e "./external/ncbi-blast+/bin/blastx" ] ||
     [ ! -e "./external/ncbi-blast+/bin/blastp" ] ); then
	packagesNotInstalled+=("BLAST+")
fi

# Check Diamond installation
if [ ! -e "./external/diamond/diamond" ]; then
	packagesNotInstalled+=("Diamond")
fi

# Check Corset installation
if [ ! -e "./external/corset/corset" ]; then
	packagesNotInstalled+=("Corset")
fi

# Check Salmon installation
if [ ! -e "./external/salmon/bin/salmon" ]; then
	packagesNotInstalled+=("Salmon")
fi

# Check TransDecoder installation
if ( [ ! -e "./external/TransDecoder/TransDecoder.LongOrfs" ] ||
     [ ! -e "./external/TransDecoder/TransDecoder.Predict" ] ); then
	packagesNotInstalled+=("TransDecoder")
fi

# Determine if any packages installed incorrectly
if [ "${#packagesNotInstalled[@]}" -eq 0 ]; then
	echo "All dependencies installed correctly."
else
	echo "ERROR: The following dependencies failed to install:"
	for package in "${packagesNotInstalled[@]}"; do
		echo "  $package"
	done
	exit 1
fi


#========================================================================
#//////////////////   BUILD SEMBLANS FROM SOURCE   \\\\\\\\\\\\\\\\\\\\\\
#========================================================================

echo "Now building Semblans ..."
mkdir ./bin/
make

# Determine if build was successful
if ( [ ! -e "./bin/preprocess" ] ||
     [ ! -e "./bin/assemble" ] ||
     [ ! -e "./bin/postprocess" ] ||
     [ ! -e "./bin/semblans" ] ); then
	echo ERROR: Semblans failed to build from source
	exit 1
fi
