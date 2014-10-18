RNAseq-Mapper
=============

To align a RNA-seq dataset towards a reference genome, it is often necessary to make use of a spliced mapper. Some popular splice mappers for RNA-seq datasets include Mapsplice, Tophat, STAR, etc. Although these mappers are generally quite fast, and works quite decently for purposes like gene expression analysis, splice junction discoveries, etc., they may not be very suitable for some niche applications like single nucleotide variant detection from RNA-seq datasets. In particular, we have noticed that these tools often lead to the identification of a significant number of false variants close to the splice sites.

To improve the accuracy of variant detection and of mapping, an alternative approach could be taken. In 



