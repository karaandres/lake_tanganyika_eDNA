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
- Separate F and R primers into 2 .csv files (format: primer name,sequence) and copy into wd
```
cp /home/kja68/lake_tanganyika/F_primers.csv ./
cp /home/kja68/lake_tanganyika/R_primers.csv ./
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
for dir in */; do mv *paired-"${dir%/}"_F.fastq.paired.fq $dir; done
for dir in */; do mv *paired-"${dir%/}"_R.fastq.paired.fq $dir; done

**Sanity check**: Check how many sequences are in each file
```
echo $(cat COI/*.fq | wc -l)/4|bc # 61092702
echo $(cat ND2/*.fq | wc -l)/4|bc # 24550130
echo $(cat *unsorted.fastq | wc -l)/4|bc # 16013784
```

**DADA2** (Callahan et al. 2016): denoise sequences with DADA2




