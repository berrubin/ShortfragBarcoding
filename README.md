# Bee barcode processing

The config.txt should contain paths to the binary exectuables for fastq-multx and USEARCH. If you want to run BLAST to obtain a rough taxonomic identification for each sample, you will also need to specify the path to a local download of NCBI's nt database. Biopython is also needed. Refer to "config.txt". 

fastq-multx can be downloaded here: https://github.com/brwnj/fastq-multx

USEARCH can be downloaded from here: http://www.drive5.com/usearch/

The i5 map file (`--i5_map`) consists of three tab-delimited columns. In the first column is the name of the plate with an individual i5 index. In the second column is the index sequence itself. And in the last column is the name of the locus being sequenced on that plate. Refer to "i5_map.txt".

The i7 map file (`--i7_map`) consists of two tab-delimited columns. In the first column is the name of the well for a particular sample and in the second column is the i7 index sequence corresponding to that well. Refer to "i7_map.txt".

Multi-threading (`--num_threads`) simply parallelizes the processing of individual samples.

For example, to demultiplex the reads in a directory named "data/", using 14 cores: 

```
shortfrag_barcoding.py -x data/read_3_index_passed_filter.fastq.gz -s data/read_2_index_passed_filter.fastq.gz \
-1 data/read_1_passed_filter.fastq.gz -2 data/read_4_passed_filter.fastq.gz -5 i5_map.txt \
-7 i7_map.txt -p 14 -c config.txt -o barcode_seqs.fa
```

If you specify an email using `-e`, then the pipeline will also attempt to assign taxonomic identifications to each sample using BLASTN. For this, you will also need to download NCBI's nucleotide database, or have a BLAST database of whatever sequences you would like to use for assigning taxonomy. The path to this database should be specified in the config.txt file with "ncbi_nt_path". If you are using a custom set of sequences then you will need to create a BLAST database using `makeblastdb`. NCBI's BLAST databases can be downloaded from their ftp site (ftp://ftp.ncbi.nlm.nih.gov/blast/db/). The procedure for this is to find the best matching sequence in the given database and then pull the taxonomy of that sequence from NCBI's server. If you're using your own database, it will need to be composed of sequences present in the NCBI database. 

The output file is just a fasta file of a sequence for each sample, assuming that a single sequence was clearly the most representative of that sample. When the identification of the most abundant sequence originating from a sample is equivocal, several sequences may be printed. If a sequence is at least half as abundant as the most abundant sequence derived from a sample, that sequence is included in the output. The number of reads that were an exact match to the given sequence is printed in the fasta header. If the BLAST taxonomy identification is used then the putative taxonomic identification for that sequence is also given in the header, along with the gi number of the closest hit sequence. If multiple sequences are given for a single sample, they are numbered consecutively.

So, for example, a fasta header may look like:  
`>MER013-A2_2;size=177;334884163;Lasioglossum_obscurum;99.35;309;313`

In this case, the sample name is "MER013-A2".  "_2" indicates that this is the second most abundant sequence from this sample that is at least half as abundant as the most abundant sequence. "size=177" means that there were 177 reads in the sample that exactly matched this sequence. "334884163" is the gi number for the best blast hit, "Lasioglossum_obscurum" is the species from which that gi number originated, "95.35" is the percent identity of the barcode sequence to the gi sequence over the length of the best HSP, which is "309" bases long. The whole barcode sequences is "313" bases long. When the best blast hit is not a Hymenopteran, this is noted in an additional field as "not_hymenoptera".

If you have another piece of information that you would like to add to the fasta file (i.e. the code for the collection site), you can create a three column tab-delimited text file and feed that into the process with `-g`. The first column should be the collections site, the second column should be the plate ID, and the third column should be the well ID. `sample_sites.txt` is an example file. If you want a header line in your file, it needs to start with a "#".