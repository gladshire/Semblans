# Semblans
Author: Miles Woodcock-Girard
*Walker Lab, UIC Department of Biological Sciences*

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; As-***Sembl***-y pipeline for tr-***ans***-criptomes

Semblans is a tool that enables the automatic assembly of *de novo* transcriptomes for non-model organisms.

Through the collation of several external packages and the leveraging of C++ data streaming performance, Semblans streamlines the necessary pre-processing, quality control, assembly, and post-assembly steps, allowing a hands-off assembly process without loss to versatility. The following diagram shows a graphical workflow of the pipeline:

![](https://live.staticflickr.com/65535/53545551413_5bd8abc933_k.jpg)

All documentation for Semblans can be found in the [wiki](https://github.com/gladshire/Semblans/wiki)

# Installation

If the user has downloaded a pre-compiled Semblans release, this step is **NOT NECESSARY**. To install Semblans, navigate to its root directory and then call:
```
./install.sh
```
Please allow several minutes for Semblans to set up the necessary packages.

By default, Semblans will not retrieve the PantherDB functional protein database for sequence annotation. **If the user intends to utilize Semblans' annotation functionality, they should instead call the following installation command:
```
./install.sh --with-panther
```
**Be aware that the PantherDB database is large (~17GB compressed; ~80GB uncompressed), and can take up to an hour for download.**

# Quick Start / Test data

Included with Semblans is a directory called 'examples'. This directory contains a very small short read dataset ("ChloroSubSet") for testing/verifying functionality of the Semblans pipeline. To test, uncompress the data from **ChloroSubSet.tar.gz**. The user should then **ensure they have a reference proteome**, as one is necessary for several of the pipeline's postprocessing stages. Links to broad, kingdom-level reference proteomes are hosted at the bottom of this document. In this example, I use the kingdom-level plant proteome. Once prepared, the user may call:
```
semblans --left ChloroSubSet_1.fq --right ChloroSubSet_2.fq --ref-proteome ensembl_plant.pep.all.fa --threads 8 --ram 10
```

Reference peptide sets for **postprocess**

[**Ensembl animal reference**](https://uofi.app.box.com/s/0rlp6q0u5uvc161mzbdr3c0xoiti63sk)

[**Ensembl plant reference**](https://uofi.app.box.com/s/lvg7x2qrxvg8hfcmue9xv4y9t1rgfb48)

[**Ensembl fungi reference**](https://uofi.box.com/s/qc4nmun4apb0pik5fqxukn4qvn3wm943)
