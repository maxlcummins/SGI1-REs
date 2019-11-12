# Salmonella genomic island 1 is widely disseminated in the Proteobacteriaceae: Identification of a potential ancestral SGI2 variant in Vibrio cholerae
A Github repository containing supplementary methodology and materials for a manuscript. DOI to be included on publication.

## SGI1-related element (SGI1-RE) dendrogram
A non-redundant nucleotide database was created by concatenating the CDS fasta files of representative SGI1-REs from the NCBI nucleotide database, including accession numbers.

This was performed using CD-HIT version 4.8.1 with the following command:

```cd-hit -i References_CDhit/SGI1-RE_cat.fasta -o CD_hit_SGI1s_100_2.fa -c 1```
