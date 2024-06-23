# Semblans
Author: Miles Woodcock-Girard
*Walker Lab, UIC Department of Biological Sciences*

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; As-***Sembl***-y pipeline for tr-***ans***-criptomes

Semblans is a tool that enables the automatic assembly of *de novo* transcriptomes for non-model organisms.

Through the collation of several external packages and the leveraging of C++ data streaming performance, Semblans streamlines the necessary pre-processing, quality control, assembly, and post-assembly steps, allowing a hands-off assembly process without loss to versatility. The following diagram shows a graphical workflow of the pipeline:

![](https://live.staticflickr.com/65535/53545551413_5bd8abc933_k.jpg)

All documentation for Semblans can be found in the [wiki](https://github.com/gladshire/Semblans/wiki)

# Installation

If the user has downloaded a pre-compiled Semblans release, this step is not necessary. To install Semblans, navigate to its root directory and then call:
```
./install.sh
```
Please allow several minutes for Semblans to set up the necessary packages.

By default, Semblans will not retrieve the PantherDB functional protein database for sequence annotation. **If the user intends to utilize Semblans' annotation functionality, they should instead call the following installation command:
```
./install.sh --with-panther
```
**Be aware that the PantherDB database is large (~17GB compressed; ~80GB uncompressed), and can take some time to download.**

# Quick Start / Test data

Included with Semblans is a directory called 'examples'. This directory contains a very small short read dataset ("ChloroSubSet") for testing/verifying functionality of the Semblans pipeline. To test, uncompress the data from **ChloroSubSet.tar.gz**. The user should then **ensure they have a reference proteome**, as one is necessary for several of the pipeline's postprocessing stages. Links to broad, kingdom-level reference proteomes are hosted at the bottom of this document. In this example, I use the kingdom-level plant proteome. Once prepared, the user may call:
```
semblans --left ChloroSubSet_1.fq --right ChloroSubSet_2.fq --prefix ChloroSubSet --ref-proteome ensembl_plant.pep.all.fa --threads 8 --ram 10
```
Some users may experience issues, particularly during the transcript assembly phase during Trinity. Common errors and solutions are hosted on our GitHub's [wiki page](https://github.com/gladshire/Semblans/wiki/Common-Errors-&-Solutions#issues-at-trinity-stage). As cataloguing these is an ongoing process, we urge users to post an issue on the Semblans repository page detailing their problem if it persists or is otherwise unaddressed by this page.

Reference peptide sets for **postprocess**

[**Ensembl animal reference**](https://www.dropbox.com/scl/fi/n49jm9i1yscrfrsj1dnq8/ensembl_animals.pep.all.fa.gz?rlkey=yemush6bm36wr4fu7dpe8h5e0&st=98kgb83l&dl=1)

[**Ensembl plant reference**](https://www.dropbox.com/scl/fi/hbvtnd9wsiwt8k7gakcmq/ensembl_plant.pep.all.fa.gz?rlkey=8cp9sn5wrt9uu4vc2pmg8xvin&st=1gesuqq0&dl=1)

[**Ensembl fungi reference**](https://www.dropbox.com/scl/fi/8as6tci331utcrqrl7txz/ensembl_fungi.pep.all.fa.gz?rlkey=eyhsv35lelnao5xgd7s9dy51e&st=oc9dz9s6&dl=1)
