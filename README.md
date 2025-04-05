# Semblans
Author: Miles Woodcock-Girard
*Walker Lab, UIC Department of Biological Sciences*

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; As-***Sembl***-y pipeline for tr-***ans***-criptomes

Semblans is a tool that enables the automatic assembly of *de novo* transcriptomes for non-model organisms. Read our application note in *Bioinformatics*.

https://doi.org/10.1093/bioinformatics/btaf003

## The easiest way to install Semblans is via Docker, or to download the latest the pre-built binaries [here](https://github.com/gladshire/Semblans/releases)

Through the integration of several external packages and the leveraging of C++ data streaming performance, Semblans streamlines the necessary pre-processing, quality control, assembly, and post-assembly steps, allowing a hands-off assembly process without loss to versatility. The following diagram shows a graphical workflow of the pipeline. The reference proteome has been omitted for simplicity, but is utilized by Diamond during the BLASTX / BLASTP steps of postprocess:

![](https://live.staticflickr.com/65535/54238413023_f9215a0ee1_o.jpg)

All documentation for Semblans can be found in the [wiki](https://github.com/gladshire/Semblans/wiki)

# Dependencies

Semblans will install most of the dependencies it requires, but make sure you have working installations of:
- [**bowtie2**](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [**jellyfish**](https://genome.umd.edu/jellyfish.html)
- [**salmon**](https://salmon.readthedocs.io/en/latest/salmon.html)
- [**samtools**](http://www.htslib.org)
- [**numpy**](https://numpy.org)

On Ubuntu this can be done by running:
```bash
sudo apt update
sudo apt install bowtie2 jellyfish salmon samtools python3-numpy
```

# Installation

The easiest way to install Semblans is via Docker, or to download the latest the binaries [here](https://github.com/gladshire/Semblans/releases).

If instead the user wishes to build from source, they must clone this repository, navigate to the Semblans root directory and then call:

```bash
./install.sh
```
Please allow several minutes for Semblans to set up the necessary packages.

By default, Semblans will not retrieve the PantherDB functional protein database for sequence annotation. **If the user intends to utilize Semblans' annotation functionality, they should instead call the following installation command:
```bash
./install.sh --with-panther
```
**Be aware that the PantherDB database is large (~17GB compressed; ~80GB uncompressed), and can take some time to download.**

# Quick Start / Test data

Included with Semblans is a directory called 'examples'. This directory contains a very small short read dataset ("ChloroSubSet") for testing/verifying functionality of the Semblans pipeline. To test, uncompress the data from **ChloroSubSet.tar.gz**. The user should then **ensure they have a reference proteome**, as one is necessary for several of the pipeline's postprocessing stages. Links to broad, kingdom-level reference proteomes are hosted at the bottom of this document. In this example, I use the kingdom-level plant proteome. Once prepared, the user may call:
```bash
semblans \
--left ChloroSubSet_1.fq \
--right ChloroSubSet_2.fq \
--prefix ChloroSubSet \
--ref-proteome ensembl_plant.pep.all.fa \
--threads 4 \
--ram 10
```
Some users may experience issues, particularly during the transcript assembly phase. Common errors and solutions are hosted on our GitHub's [wiki page](https://github.com/gladshire/Semblans/wiki/Common-Errors-&-Solutions#issues-at-trinity-stage). As cataloguing these is an ongoing process, we urge users to post an issue on the Semblans repository page detailing their problem if it persists or is otherwise unaddressed by this page.

Reference peptide sets (gzipped FASTA) for the **postprocess** step:

**Ensembl animal reference** (3.1 GiB) [[Option 1](https://uofi.box.com/shared/static/0rlp6q0u5uvc161mzbdr3c0xoiti63sk) | [Option 2](https://www.dropbox.com/scl/fi/n49jm9i1yscrfrsj1dnq8/ensembl_animals.pep.all.fa.gz?rlkey=yemush6bm36wr4fu7dpe8h5e0&st=98kgb83l&dl=1)]

**Ensembl fungi reference** (4.7 GiB) [[Option 1](https://uofi.box.com/shared/static/qc4nmun4apb0pik5fqxukn4qvn3wm943) | [Option 2](https://www.dropbox.com/scl/fi/8as6tci331utcrqrl7txz/ensembl_fungi.pep.all.fa.gz?rlkey=eyhsv35lelnao5xgd7s9dy51e&st=oc9dz9s6&dl=1)]

**Ensembl plant reference** (2.1 GiB) [[Option 1](https://uofi.box.com/shared/static/lvg7x2qrxvg8hfcmue9xv4y9t1rgfb48) | [Option 2](https://www.dropbox.com/scl/fi/hbvtnd9wsiwt8k7gakcmq/ensembl_plant.pep.all.fa.gz?rlkey=8cp9sn5wrt9uu4vc2pmg8xvin&st=1gesuqq0&dl=1)]

# Citing

Woodcock-Girard, M. D., Bretz, E. C., Robertson, H. M., Ramanauskas, K., Hampton-Marcell, J. T., & Walker, J. F. (2024). Semblans: automated assembly and processing of RNA-seq data. Bioinformatics (Oxford, England), 41(1), btaf003. https://doi.org/10.1093/bioinformatics/btaf003
