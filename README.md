# aDNA_GenoCaller
aDNA_GenoCaller is a python program that calls genotypes from bam files at positions/regions specified in a bed file while taking into account post mortem damage as estimate by MapDamage

This version of the program will take any genotype with a low quality heterozygote call (Q<30) and convert to the next best homozygote call. pysam and matplotlib need to be installed to run the script. 

Three files are created:
An emit all vcf file noting the call for all base pairs in the bed file.
A vcf file with only those sites with evidence for at least one alternative allele and that passes a predetermined QUAL filter.
A haploid emit all vcf (ish), that gives you the most likely base under a haploid model. If two or more basepairs are tied, the reported allele is randomly chosen.

To run aDNA_GenoCaller type:

./aDNA_GenoCaller.py indexed_bamfile bed_file reference_genome 5CT_mapdamage_file 3GA_mapdamage_file
