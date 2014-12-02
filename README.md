RNAseq-Mapper
=============

To align a RNA-seq dataset towards a reference genome, it is often necessary to make use of a spliced mapper. Some popular splice mappers for RNA-seq datasets include Mapsplice, Tophat, STAR, etc. Although these mappers are generally quite fast, and works quite decently for purposes like gene expression analysis, splice junction discoveries, etc., they may not be very suitable for some niche applications like single nucleotide variant detection from RNA-seq datasets. In particular, we have noticed that these tools often lead to the identification of a significant number of false variants close to splice sites.

To improve the accuracy of variant detection and of mapping, we took an alternative approach by using a highly accurate mapper (BWA/BWA-MEM) to map the RNA-seq datasets to a human genome in combination with a spliced junction database. Outline here is the overall methodology for performing such an alignment together with some of my code for the entire pipeline.

## Overview
This is the code for the RNA-seq aligner using a hg19 reference genome
as well as a junction database of sequences


## Pipeline overview
The pipeline can be divided into a few parts

1. Generation of a junction database sequences based on known annotations (UCSC, ENSEMBL, GENECODE, Refseq).

2. Alignment of RNA-seq data to hg19+junction_database reference fasta file concurrently.

3. Conversion of junction mapped reads onto hg19 coordinates.

4. Selection of best pair among the reads.

(Only code for Part 1 is provided. Please contact me if you need the code of Parts 3 and 4.)
