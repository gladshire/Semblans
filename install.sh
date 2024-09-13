#!/usr/bin/env bash
####################################################
# Suppress printing of error messages
# exec 2>/dev/null
# Stop on first error
set -o errexit
# Set trap on ERR to be inherited by shell functions
set -o errtrace
# Trap errors
trap 'echo Error at line: $LINENO' ERR
####################################################

# Determine whether to install PANTHER for annotation
install_panther=false
for flag in "$@"; do
  if [ "$flag" == "--with-panther" ] || [ "$flag" == "-p" ]; then
    install_panther=true
  fi
  if [ "$flag" == "--help" ]; then
    echo "  ./install.sh [--with-panther / -p]"
    exit 0
  fi
done

echo "Initiating install of Semblans"

if $install_panther
then
  echo "PANTHER will be installed. Please allow several extra minutes for its download."
  echo "If user does not wish to perform annotations, they may omit the --with-panther / -p flag."
else
  echo "PANTHER will NOT be installed, and is required for transcript annotation."
  echo "If user wishes to perform annotations, they should include the --with-panther / -p flag."
fi

# Prepare Semblans directory structure
rm -rf ./include ./lib ./external ./data ./bin

mkdir -p ./include
mkdir -p ./lib
mkdir -p ./external

if $install_panther; then
	mkdir -p ./data
fi

#========================================================================
#/////////////   INSTALLATION OF THIRD-PARTY LIBRARIES   \\\\\\\\\\\\\\\\
#========================================================================

echo "Now installing required libraries"

#
# ToDo: Miles, boost failed to build for me because of this:
#       "--with-python=python3"    <-- with Python 3.12 as my default
#                                      Python install.
#       "--with-python=python3.10" <-- worked fine.
#
# Install boost libraries
if  [ ! -f ./lib/libboost_filesystem.a ] ||
    [ ! -f ./lib/libboost_regex.a ] ||
    [ ! -f ./lib/libboost_system.a ] ||
    [ ! -f ./lib/libboost_locale.a ]; then
	echo "  Installing Boost libraries ..."
	wget -q https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
	tar -xf boost_1_81_0.tar.gz
	cd boost_1_81_0 || return 1
	{
    ./bootstrap.sh --prefix=../ --with-python=python3.10
	./b2 install cxxflags="-std=c++11" link=static
    } ||
    {
    ./bootstrap.sh --prefix=../ --with-python=python3
    ./b2 install cxxflags="-std=c++11" link=static
    }
	mv LICENSE_1_0.txt ../include/boost/
	cd ..
	rm -rf boost_1_81_0*
else
	echo "Boost already installed. Skipping ..."
fi

# Install rapidXML
if [ ! -f ./include/rapidxml/rapidxml.hpp ]; then
	echo "  Installing rapidXML library ..."
	wget -q https://sourceforge.net/projects/rapidxml/files/rapidxml/rapidxml%201.13/rapidxml-1.13.zip
	unzip rapidxml-1.13.zip -d ./include/
	mv ./include/rapidxml-1.13 ./include/rapidxml
	rm rapidxml-1.13.zip
else
	echo "RapidXML already installed. Skipping ..."
fi

# Install libconfini
echo "  Installing libconfini library ..."
wget -q https://github.com/madmurphy/libconfini/releases/download/1.16.4/libconfini-1.16.4-x86_64-bin.tar.xz
tar -xf libconfini-1.16.4-x86_64-bin.tar.xz
mkdir -p ./include/libconfini
mv ./usr/include/* ./include/libconfini/
mv ./usr/lib/* ./lib/
mkdir -p ./lib/pkgconfig
mv ./usr/lib/pkgconfig/* ./lib/pkgconfig/
mv ./usr/share/doc/libconfini/AUTHORS ./include/libconfini/
mv ./usr/share/doc/libconfini/COPYING ./include/libconfini/
rm libconfini-1.16.4-x86_64-bin.tar.xz
rm -rf ./usr

# Install libcurl
echo " Installing libcurl library ..."
wget -q https://curl.se/download/curl-8.1.2.tar.gz
tar -xf curl-8.1.2.tar.gz
mkdir -p ./include/curl
mv ./curl-8.1.2/include/curl/* ./include/curl
mv ./curl-8.1.2/COPYING ./include/curl
rm -rf curl-8.1.2
rm curl-8.1.2.tar.gz

#========================================================================
#/////////////   INSTALLATION OF THIRD-PARTY PACKAGES    \\\\\\\\\\\\\\\\
#========================================================================

echo "Now installing required packages ..."

# Install HMMER
if [ ! -e "./external/hmmer/bin/hmmscan" ]; then
	echo "Installing HMMER ..."
	wget -q http://eddylab.org/software/hmmer/hmmer.tar.gz
	tar -xf hmmer.tar.gz -C ./external/
	mv ./external/hmmer* ./external/hmmer
	cd ./external/hmmer || return 1
	./configure --prefix "$PWD"
	make
	make install
	cd ../../
	rm hmmer.tar.gz
else
	echo "HMMER already installed. Skipping ..."
fi

# Install PANTHER HMM Scoring components
if $install_panther ; then
	if [ ! -e "./external/pantherScore/pantherScore2.2.pl" ]; then
  		echo "Downloading PANTHER components ..."
		wget -q http://data.pantherdb.org/ftp/hmm_scoring/current_release/PANTHER18.0_hmmscoring.tgz
		tar -xzf PANTHER*.tgz
		mv target panther_db/
		mv panther_db ./data/
		rm PANTHER*.tgz
		wget -r -q --no-parent http://data.pantherdb.org/ftp/hmm_scoring/17.0/pantherScore2.2/
		mv data.pantherdb.org/ftp/hmm_scoring/17.0/pantherScore2.2 ./external/pantherScore
		chmod +x ./external/pantherScore/pantherScore2.2.pl
		rm -rf data.pantherdb.org
	else
		echo "PANTHER components found. Skipping ..."
	fi
fi


# Install NCBI sra-tools
if  [ ! -e "./external/sratoolkit/bin/prefetch" ] &&
    [ ! -e "./external/sratoolkit/bin/fasterq-dump" ]; then
	echo "Installing SRA-Tools ..."
	wget -q --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
	tar -xf sratoolkit.tar.gz -C ./external/
	mv ./external/sratoolkit* ./external/sratoolkit
	cd ./external/sratoolkit || return 1
	wget -q https://raw.githubusercontent.com/ncbi/sra-tools/master/LICENSE
	cd ../../
	rm sratoolkit.tar.gz
else
	echo "SRA-Tools already installed. Skipping ..."
fi

# Install FastQC
if [ ! -e "./external/FastQC/fastqc" ]; then
	echo "Installing FastQC ..."
	wget -q --output-document fastqc.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
	unzip fastqc.zip -d ./external/
	chmod 755 ./external/FastQC/fastqc
	rm fastqc.zip
else
	echo "FastQC already installed. Skipping ..."
fi

# Install Rcorrector
if [ ! -e "./external/Rcorrector/rcorrector" ]; then
	echo "Installing Rcorrector ..."
	wget -q --output-document rcorrector.tar.gz https://github.com/mourisl/Rcorrector/archive/refs/tags/v1.0.7.tar.gz
	tar -xf rcorrector.tar.gz -C ./external/
	mv ./external/Rcorrector* ./external/Rcorrector
	cd ./external/Rcorrector || return 1
	make
	cd ../..
	rm rcorrector.tar.gz
else
	echo "Rcorrector already installed. Skipping ..."
fi

# Install Trimmomatic
if [ ! -e "./external/Trimmomatic/trimmomatic-0.39.jar" ]; then
	echo "Installing Trimmomatic ..."
	wget -q --output-document trimmomatic.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
	unzip trimmomatic.zip -d ./external/
        mv ./external/Trimmomatic* ./external/Trimmomatic
	cd ./external/Trimmomatic/adapters || return 1
	echo -e "$(cat NexteraPE-PE.fa)\n$(cat TruSeq2-PE.fa)\n$(cat TruSeq2-SE.fa)\n$(cat TruSeq3-PE-2.fa)\n$(cat TruSeq3-PE.fa)\n$(cat TruSeq3-SE.fa)" > TruSeq_all.fa
	cd ../../..
	rm trimmomatic.zip
else
	echo "Trimmomatic already installed. Skipping ..."
fi

# Install Kraken2
if [ ! -e "./external/kraken2/kraken2" ]; then
	echo "Installing Kraken2 ..."
	wget -q --output-document kraken2.tar.gz https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.tar.gz
	tar -xf kraken2.tar.gz -C ./external/
	mv ./external/kraken2-2.1.2 ./external/kraken2
	cd ./external/kraken2 || return 1
	./install_kraken2.sh .
	cd ../..
	rm kraken2.tar.gz
else
	echo "Kraken2 already installed. Skipping ..."
fi

# Install Trinity
if [ ! -e "./external/trinityrnaseq/Trinity" ]; then
	echo "Installing Trinity ..."
	wget -q --output-document trinity.tar.gz https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz
	tar -xf trinity.tar.gz -C ./external/
	mv ./external/trinityrnaseq-v2.15.1 ./external/trinityrnaseq
	cd ./external/trinityrnaseq/trinity-plugins/bamsifter || return 1
	sed -i '1s/^/#include <string> \n/' sift_bam_max_cov.cpp
	cd ../..
	make
	cd ../..
	rm trinity.tar.gz
else
	echo "Trinity already installed. Skipping ..."
fi

# Install NCBI BLAST
if [ ! -e "./external/ncbi-blast+/bin/blastx" ] ||
   [ ! -e "./external/ncbi-blast+/bin/blastp" ]; then
	echo "Installing NCBI-BLAST+ ..."
	wget -q --output-document ncbi-blast.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz
	tar -xf ncbi-blast.tar.gz -C ./external/
	mv ./external/ncbi-blast* ./external/ncbi-blast+
	rm ncbi-blast.tar.gz
else
	echo "NCBI-BLAST+ already installed. Skipping ..."
fi

# Install Diamond
if [ ! -e "./external/diamond/diamond" ]; then
	echo "Installing Diamond ..."
	mkdir -p ./external/diamond-2.1.7
	wget -q --output-document diamond.tar.gz https://github.com/bbuchfink/diamond/releases/download/v2.1.7/diamond-linux64.tar.gz
	tar -xf diamond.tar.gz -C ./external/diamond-2.1.7/
	mv ./external/diamond-2.1.7 ./external/diamond
	cd ./external/diamond || return 1
	wget -q https://raw.githubusercontent.com/bbuchfink/diamond/master/LICENSE
	cd ../..
	rm diamond.tar.gz
else
	echo "Diamond already installed. Skipping ..."
fi


# Install Corset
if [ ! -e "./external/corset/corset" ]; then
	echo "Installing Corset ..."
	wget -q --output-document corset.tar.gz https://github.com/Oshlack/Corset/releases/download/version-1.09/corset-1.09-linux64.tar.gz
	tar -xf corset.tar.gz -C ./external/
	mv ./external/corset-1.09-linux64 ./external/corset
	cd ./external/corset || return 1
	chmod -x LICENSE
	chmod -x COPYING
	cd ../..
	rm corset.tar.gz
else
	echo "Corset already installed. Skipping ..."
fi

# Install STAR
# if [ ! -e "./external/STAR/bin" ]; then
# 	echo "Installing STAR ..."
# 	wget -q --output-document star.tar.gz https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11a.tar.gz
# 	tar -xf star.tar.gz -C ./external/
# 	mv ./external/STAR* ./external/STAR
# 	cd ./external/STAR/source || return 1
# 	make STAR
# 	cd ../../..
# 	rm star.tar.gz
# else
# 	echo "STAR alread installed. Skipping ..."
# fi

# Install Salmon
if [ ! -e "./external/salmon/bin/salmon" ]; then
	echo "Installing Salmon ..."
	wget -q --output-document salmon.tar.gz https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
	tar -xf salmon.tar.gz -C ./external/
	mv ./external/salmon* ./external/salmon
	cd ./external/salmon || return 1
	wget -q https://raw.githubusercontent.com/COMBINE-lab/salmon/master/LICENSE
	cd ../..
	rm salmon.tar.gz
else
	echo "Salmon already installed. Skipping ..."
fi

# Install TransDecoder
if [ ! -e "./external/TransDecoder/TransDecoder.LongOrfs" ] ||
   [ ! -e "./external/TransDecoder/TransDecoder.Predict" ]; then
	echo "Installing TransDecoder ..."
	wget -q --output-document transdecoder.tar.gz https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.0.tar.gz
	tar -xf transdecoder.tar.gz -C ./external/
	mv ./external/TransDecoder-TransDecoder-v5.7.0 ./external/TransDecoder
	rm transdecoder.tar.gz
else
	echo "TransDecoder already installed. Skipping ..."
fi

#========================================================================
#//////////   TEST SUCCESSFUL INSTALLATION OF DEPENDENCIES    \\\\\\\\\\\
#========================================================================

packagesNotInstalled=()

# Check HMMER installation
if [ ! -e "./external/hmmer/bin/hmmscan" ]; then
	packagesNotInstalled+=("HMMER")
fi

# Check pantherScore installation
if $install_panther ; then
  if [ ! -e "./external/pantherScore/pantherScore2.2.pl" ]; then
    packagesNotInstalled+=("pantherScore")
  fi
fi

# Check SRA Toolkit installation
if [ ! -e "./external/sratoolkit/bin/prefetch" ] &&
   [ -e "./external/sratoolkit/bin/fasterq-dump" ]; then
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
if [ ! -e "./external/ncbi-blast+/bin/blastx" ] ||
   [ ! -e "./external/ncbi-blast+/bin/blastp" ]; then
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
if [ ! -e "./external/TransDecoder/TransDecoder.LongOrfs" ] ||
   [ ! -e "./external/TransDecoder/TransDecoder.Predict" ]; then
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
mkdir -p ./bin
make

# Determine if build was successful
if [ ! -e "./bin/preprocess" ] ||
   [ ! -e "./bin/assemble" ] ||
   [ ! -e "./bin/postprocess" ] ||
   [ ! -e "./bin/semblans" ]; then
	echo "ERROR: Semblans failed to build from source"
	exit 1
else
	echo "Semblans built successfully"
fi
