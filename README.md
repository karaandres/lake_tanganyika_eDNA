# Lake Tanganyika eDNA metabarcoding
## This repository contains the files related to the data processing and data analysis for eDNA metabarcoding of fishes in Lake Tanganyika. 

### Overview
Add project overview here (sampling info, library prep, sequencing details).
* COI primers
> COI_F:	GGTACTGGTTGAACTGTTTATCCTCC

> COI_R:	TANACYTCNGGRTGNCCRAARAAYCA

* ND2 primers
> ND2_F:	CATACCCCAARMATGTTGGTTAAA

> ND2_R:	GCAAGYTTTTGTCAGGTTGARAG

**Note**: COI and ND2 libraries were sequenced together, where each sample had a unique N5/N7 barcode and sequences were parsed into 2 loci during bioinformatic processing (see below). 

### Bioinformatic processing steps
- Create directory on HPC, move all sequence files into directory, clone this repo
```
mkdir /workdir/kja68/ && cd /workdir/kja68/
cp /home/kja68/lake_tanganyika/raw_sequence_files/*.gz ./
git clone https://github.com/karaandres/lake_tanganyika_eDNA
```
**Trimmomatic**: Trim reads that match Nextera adapters or are shorter than 35 bp
- Input is raw reads demultiplexed by sample (.fastq.gz), output is trimmed reads (paired.fastq.gz and unpaired.fastq.gz)
- Make script executable, run on all files
```
cd lake_tanganyika_eDNA/scripts/
chmod u+x trimmomatic_loop.sh
cd /workdir/kja68/
lake_tanganyika_eDNA/scripts/trimmomatic_loop.sh
```

**Split_on_Primer**: Split reads into separate files by primer
- Input is paired trimmed reads (paired.fastq.gz), output is separate files for each locus (COI and ND2)
- Clone script, make executable
```
git clone https://github.com/marcomeola/Split_on_Primer.git
chmod u+x Split_on_Primer/src/Split_on_Primer_fixed.py
```
- Copy F and R primer files into directory (format: primer name,sequence)
```
cp lake_tanganyika_eDNA/files/F_primers.csv ./
cp lake_tanganyika_eDNA/files/R_primers.csv ./
```
- Make directory for split sequence files and move trimmed files
```
mkdir /workdir/kja68/split_files/
cp *_paired.fastq.gz /workdir/kja68/split_files/
```
- Unzip files (code only works on fasta or fastq) 
```
gzip -d /workdir/kja68/split_files/*_paired.fastq.gz
```
- Run script separately for F and R reads; separated files will appear in `split_files/` directory
- Separated files will get an extension on their file name with the locus name, e.g. filename-locusname_F.fastq.gz
- Might take a while for several loci -- use screen so the session stays active
```
screen
for file in split_files/*_R1_*; do Split_on_Primer/src/Split_on_Primer_fixed.py -f "$file" -p F_primers.csv -m 2; done
for file in split_files/*_R2_*; do Split_on_Primer/src/Split_on_Primer_fixed.py -f "$file" -p R_primers.csv -m 2; done
```

**fastq-pair**: Previous script did not provide option to retain matching F and R reads in identical order (required for downstream analyses)
- Input is trimmed files demultiplexed by locus (locusname.fastq.gz), output is matched F and R reads (paired.fq)
- Clone script to match forward and reverse reads and make executable
```
git clone https://github.com/linsalrob/fastq-pair
cd fastq-pair/
mkdir build && cd build
gcc -std=gnu99   ../main.c ../robstr.c ../fastq_pair.c ../is_gzipped.c  -o fastq_pair
```
- Change directory and run script
```
cd split_files/
for f1 in *_F.fastq; do f2a=${f1%%_F.fastq}"_R.fastq"; f2=${f2a/R1-001.fastq.gz/R2-001.fastq.gz}; ../fastq-pair/build/fastq_pair -t 10000 $f1 $f2; done
```
- Make directories for each locus (use the same locus names that are appended to your file names) and move split files into the appropriate directory by locus name 
```
mkdir COI ND2
```
- Loops thru the locus directories you just created and grab files by locus name
```
for dir in */; do mv *paired-"${dir%/}"_F.fastq.paired.fq $dir; done
for dir in */; do mv *paired-"${dir%/}"_R.fastq.paired.fq $dir; done
```

**Sanity check**: Check how many sequences are in each file
```
echo $(cat COI/*.fq | wc -l)/4|bc # 61092702
echo $(cat ND2/*.fq | wc -l)/4|bc # 24550130
echo $(cat *unsorted.fastq | wc -l)/4|bc # 16013784
```

### Denoise sequences with DADA2 (Callahan et al. 2016): 
- Load packages
```
/programs//R-4.0.5/bin/R
library(ShortRead); packageVersion("ShortRead") # 1.48.0
library(dada2); packageVersion("dada2") # 1.19.2
```
- Run separately for COI and ND2
- Input is F and R sequences demultiplexed by sample and locus and matched in order
- Output is ASV by sample matrix (number of reads per ASV in each sample)
- Some extra steps added because we have lots of samples with no reads (blanks) and DADA2 doesn't like that
```
  path <- paste("/workdir/kja68/lake_tanganyika/split_files/COI", loci[2],sep="")
  setwd(path)
  fnFs <- sort(list.files(path, pattern="F.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="R.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_F.fastq.gz
	sample.names <- sapply(strsplit(basename(fnFs), "-R1"), `[`, 1)

# Inspect read quality profiles -- only for a subset of samples (1:20) and for those with lengths > 0
	fnFs_qual <- NULL
	fnRs_qual <- NULL
	for (i in 1:20){
	srq <- readFastq(fnFs[i])
	seqlen.tab <- table(width(srq))
	if (nrow(seqlen.tab)>0)
	fnFs_qual <- c(fnFs_qual,i)
	srq <- readFastq(fnRs[i])
	seqlen.tab <- table(width(srq))
	if (nrow(seqlen.tab)>0)
	fnRs_qual <- c(fnRs_qual,i)
	}

	pdf("QualityProfile.pdf")
	plotQualityProfile(fnFs[fnFs_qual])
	plotQualityProfile(fnRs[fnRs_qual])
	dev.off()

# Place filtered files in filtered/ subdirectory
	filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
	filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
	names(filtFs) <- sample.names
	names(filtRs) <- sample.names

# Filter and trim
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,290), trimLeft = c(26,26),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE, verbose=TRUE)
	head(out)

# Check how many samples passed filter and remove files that did not pass the filter
	table(file.exists(filtFs))
	table(file.exists(filtRs))
	filtFs <- filtFs[file.exists(filtFs)]
	filtRs <- filtRs[file.exists(filtRs)]

# Learn errors
	errF <- learnErrors(filtFs, multithread=TRUE)
	errR <- learnErrors(filtRs, multithread=TRUE)
	pdf("Error_plot.pdf")
	plotErrors(errF, nominalQ=TRUE)
	dev.off()
	dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
	dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
	mergers_20_1 <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap=20, maxMismatch=1, verbose=TRUE)

# Construct ASV table 
	seqtab <- makeSequenceTable(mergers_20_1)
	dim(seqtab)

# Inspect distribution of sequence lengths (ASV table)
	table(nchar(getSequences(seqtab)))
	write.csv(table(nchar(getSequences(seqtab))), "Sequence_lengths.csv")

# remove chimeras 
	seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
	sum(seqtab.nochim)/sum(seqtab)
	dim(seqtab.nochim)

# Save the tables and workspace
	write.csv(t(seqtab), "seqtab.csv")
	write.csv(t(seqtab.nochim), "seqtab.nochim.csv")
	save.image(file='Dada2.RData')

# Save fasta files
	uniquesToFasta(getUniques(seqtab), fout="uniqueSeqs.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab)))))
	uniquesToFasta(getUniques(seqtab.nochim), fout="uniqueSeqs.nochim.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab.nochim)))))
	head(out)
	mergers <- mergers_20_1

# Track reads through the pipeline
	getN <- function(x) sum(getUniques(x))
	track <- cbind(sum(out[,1]), sum(out[,2]), sum(sapply(dadaFs, getN)), sum(sapply(dadaRs, getN)), sum(sapply(mergers, getN)), sum(rowSums(seqtab.nochim)))
	colnames(track) <- c("input", "filtered", "denoisF", "denoisR", "merged", "nonchim")
	track

	pdf("Track_reads.pdf")
	barplot(colSums(track))
	dev.off()
```

### BLAST: Assign taxonomy to ASVs
- Copy the nucleotide and taxonomic database to your workdir
```
cp /shared_data/genome_db/BLAST_NCBI/nt* ./
cp /shared_data/genome_db/BLAST_NCBI/taxdb.* ./
```
- Run blastn, output is fmt 7 with taxonomic ID info (staxids)
```
blastn -query COI/uniqueSeqs.nochim.fasta -db nt -outfmt '7 qseqid saccver pident evalue qstart qend length sscinames staxids' -evalue 1e-3 -max_target_seqs 5 -num_threads 24 -out blstn_nonchim_fmt7_LT_COI.txt
blastn -query ND2/uniqueSeqs.nochim.fasta -db nt -outfmt '7 qseqid saccver pident evalue qstart qend length sscinames staxids' -evalue 1e-3 -max_target_seqs 5 -num_threads 24 -out blstn_nonchim_fmt7_LT_ND2.txt
```
- Download and unncompress higher taxnomic database from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
```
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -zxvf taxdump.tar.gz
```
- Match staxids to taxonomic information: use tax_trace.pl to match qseqid and staxids to taxonomic database
- Run separately for COI and ND2
```
git clone https://github.com/theo-allnutt-bioinformatics/scripts
chmod u+x scripts/tax_trace.pl
cut -f1,9 blstn_nonchim_fmt7_LT_COI.txt > taxids_COI.txt
perl scripts/tax_trace.pl taxdump/nodes.dmp taxdump/names.dmp taxids_COI.txt taxids_COI_export.txt
cut -f1,9 blstn_nonchim_fmt7_LT_ND2.txt > taxids_ND2.txt
perl scripts/tax_trace.pl taxdump/nodes.dmp taxdump/names.dmp taxids_ND2.txt taxids_ND2_export.txt
```

### Analyses (coming soon)
