#!/bin/env bash

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
if ( [ ! -e "./external/ncbi-blast+/blastx" ] ||
     [ ! -e "./external/ncbi-blast+/blastp" ] ); then
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
		echo "$package"
	done
	exit 1
fi
