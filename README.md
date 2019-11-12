# Salmonella genomic island 1 is widely disseminated in the Proteobacteriaceae: Identification of a potential ancestral SGI2 variant in Vibrio cholerae
A Github repository containing supplementary methodology and materials for a manuscript. DOI to be included on publication.

## SGI1-related element (SGI1-RE) dendrogram
A non-redundant nucleotide database was created by concatenating the CDS fasta files of representative SGI1-REs from the NCBI nucleotide database, including accession numbers.

This was performed using CD-HIT version 4.8.1 with the following command:

```cd-hit -i References_CDhit/SGI1-RE_cat.fasta -o CD_hit_SGI1s_100.fa -c 1```

Subsequently, BLASTn was used to screen samples under investigation for the carriage of genes in the non-redundant nucleotide database:

```for f in assemblies/*.fasta; do blastn -num_threads 2 -evalue 0.001 -db CD_hit_SGI1s_100.fa -query ${f} -out output/${f}.out -outfmt "6 qseqid stitle sseqid pident length slen sstart send qstart qend qlen mismatch gapopen evalue bitscore"; done```

Following this, BLAST data for each assembly was tagged with the filename and concatenated into a single file:

``` cd output;
    for f in *.out;
    do paste $f <(yes $f | head -n $(cat $f | wc -l)) > $f.new;
    done;
    cat *.new > CD-hit_SGI1-REs.txt```
