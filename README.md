# Simulate SNPs

---

## Methodology
1. Generate random SNP positions
2. Filter random positions
3. Select 10,000 SNPs
4. Randomly select ALT alleles
5. Create an ALT Fasta file
6. Simulate Illumina Paired-End Reads for REF and ALT Fastas
7. Pysim mixture: define somatic allele fractions
8. Align mixed FASTQs with BWA-MEM 
9. Subsample BAMs to desired coverage

---


## 1. Generate random SNPs

```
$ randomBed -l 1 -n 1000000 -g human.hg19.22.genome  | cut -f 1-3  | sortBed >sim_snps.22.bed
```

## 2 Remove SNPs within excluded regions

### Excluded Regions:

* RepeatMasker 
* Short Tandem Repeats (STR)
* Segmental Duplications (SegDups)
* Unmappable Regions

```
$ intersectBed -a sim_snps.22.bed -b repeatMasker_hg19.bed.gz \
                                     hg19_str.bed.gz \
                                     hg19_segdup.bed \
                                     hg19_umap_k24_mask.bed.gz \
                                     -wa -v >sim_snps.22.filtered.bed
```

## 3. Randomly select 10,000 SNPs

```
$ shuf sim_snps.22.filtered.bed | head -n 10000 | sortBed >sim_snps.22.filtered.10k.bed
```

## 4. Make a VCF with ALT alleles

```
$ python src/make_snps.py sim_snps.22.filtered.10k.bed
```

Output located in `sim_snps.22.filtered.vcf`

```
$ head sim_snps.22.filtered.vcf 

##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	17280662	SIMSNP1	T	A	.	.	.
22	17280665	SIMSNP2	A	G	.	.	.
22	17280731	SIMSNP3	G	T	.	.	.
22	17280932	SIMSNP4	C	G	.	.	.
22	17283184	SIMSNP5	A	T	.	.	.
22	17289733	SIMSNP6	A	C	.	.	.
22	17292475	SIMSNP7	C	G	.	.	.
22	17305169	SIMSNP8	A	C	.	.	.
```

`make_snps.py` takes as input a BED file and randomly selects an ALT allele. 

Line 5 of the script requires the full path of a `samtools faidx` indexed FASTA file 
```
fasta = pysam.FastaFile("/home/dantakli/ref/human_g1k_v37.fasta")
```
The script determines the reference nucleotide for a given position and randomly selects an ALT allele different from the reference. 

For example, at chr22:17280662 the reference allele is `T`. The script will randomly pick `A`,`G`, or `C` as the ALT allele.

## 5. GATK FastaAlternateReferenceMaker

Use GATK to create a FASTA file containing ALT alleles. This FASTA file will be used to simulate Illumina reads.

First generate a FASTA file of the contig of interest (here chromosome 22)

```
$ samtools faidx human_g1k_v37.fasta 22 >human_g1k_v37.22.fasta
$ samtools faidx human_g1k_v37.22.fasta # index the FASTA file
$ bwa index human_g1k_v37.22.fasta # bwa index the FASTA file for alignment steps
```

Generate the ALT FASTA file:

*Edit this script to define output file* 

```
bash scripts/gatk_fasta_alt_maker.sh in.vcf
```

**GATK command:**

```
java -Xmx16g -jar GenomeAnalysisTK.jar \
                  -T FastaAlternateReferenceMaker \
                  -R human_g1k_v37.22.fasta \
                  -o human_g1k_v37.22.ALT.fasta \
                  -L 22 -V sim_snps.22.filtered.vcf
```

## 6. Simulate Illumina Paired-End reads with ART 

This step uses the pysim script `run_art.py`

[Install this fork of Pysim (from a rotation student who fixed pysim bugs)](https://github.com/auberginekenobi/pysim)

```
git clone --recursive https://github.com/auberginekenobi/pysim.git
```

[ART Linux source:](https://www.niehs.nih.gov/research/resources/assets/docs/artsrcmountrainier2016.06.05linux.tgz)

[ART website](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
  
* REF commands

```
$ python pysim/run_art.py -i human_g1k_v37.22.fasta -l 150 -f 1000 -m 400 -s 90 -o sim_snps.REF.22.1000x
```

* ALT commands

```
$ python pysim/run_art.py -i human_g1k_v37.22.ALT.fasta -l 150 -f 1000 -m 400 -s 90 -o sim_snps.ALT.22.1000x
```

The ART commands will generate two FASTQ output files: `sim_snps.ALT.22.1000x.1.fq` and `sim_snps.ALT.22.1000x.2.fq`

ART options:

* `-l` : Read length (bp)
* `-f` : fold coverage (1000x). 
* `-m` : mean size of the DNA fragments for paired-end simulations
* `-s`: standard deviation of the RNA gragment size for paired-end simulations (I think 90 is default)

**Note on coverage**
Due to possible downsampling observed after the subsequent pysim mixture step, I recommend simulating higher coverage and down-sampling to the desired coverage after the mixture step. 

I found 1000x simulations are around 670X after the mixture step with bedtools genomecov. This might be due to ART simulating "N" nucleotides found in unmappable regions of the genome (The p-arm of chr22 for example). 

## 7. pysim mixture_v1.py

**Note on mixture_v1.py**
This script is a memory hog. I found it better to split my FASTQs into files of 100,000,000 lines first. 

```
$ split -l 100000000 -d sim_snps.ALT.22.1000x.1.fq sim_snps.ALT.22.1000x.1.
# files sim_snps.ALT.22.1000x.1.00 sim_snps.ALT.22.1000x.1.01 
# sim_snps.ALT.22.1000x.1.02 sim_snps.ALT.22.1000x.1.03 are created
```
I then use a perl script I wrote to format the file names for pysim.

```
$ perl src/rename_split_fqs.pl
```

You'll need to change the path to the files in the script. The script will print out `mv` commands that can be written to a shell script.

**Note on pysim conf format**
Pysim uses `.conf` files to mix the FASTQs. In this file you'll define the degree of mixture and the input file prefixes.

The input prefix for FASTQ file must include all text before the paired-end id <prefix.1.fq>,<prefix.2.fq>. 

Each row is a FASTQ to mix. There are two columns delimited by tabs: the first is the FASTQ prefix, the second is the degree of mixture. 

For example: 1% Allele Fraction
```
/oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/split/ref.22.1000x.00.    0.99
/oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/split/alt.22.1000x.00.    0.01
```

To generate the `.conf` file and pysim commands, you can run this script

```
$ perl src/pysim_mixture_cmd.pl
```

**Pysim mixture command**

```
$ python pysim/mixture_v1.py -i <in.conf> -o <out.FASTQ>
```

## 8. bwa mem align

This is my pipeline. I also have a perl script that will generate shell scripts to run in a job-array
```
$ perl src/bwa_snp_cmd.pl
```

**Note**: before running `bwa mem` ensure you ran `bwa index` on your FASTA file.

**bwa-mem pipeline**
```
# align 
bwa mem -t <threads> -c 500 human_g1k_v37.22.fasta snp_sim.1000x.01AF.1.fq snp_sim.1000x.01AF.2.fq -R "@RG\tID:01AF\tSM:01AF\tLB:01AF\tPL:ILLUMINA" | samtools view -bh -@ <threads> -O BAM -o snp_sim.1000x.01AF.bwamem.bam

# sort
mkdir -p snp_sim.1000x.01AF_tmp/
samtools sort -@ <threads> -o snp_sim.1000x.01AF.bwamem.sorted.bam -O BAM -T snp_sim.1000x.01AF_tmp/snp_sim.1000x.01AF snp_sim.1000x.01AF.bwamem.bam
rm -rf snp_sim.1000x.01AF_tmp/

# index
samtools index snp_sim.1000x.01AF.bwamem.sorted.bam

# check coverage
bedtools genomecov -ibam snp_sim.1000x.01AF.bwamem.sorted.bam -g human.hg19.22.genome >sim_snp.1000x.01AF.genomecov.txt
```


## 9. Sub-sample BAM

You can calculate the coverage `bedtools genomecov` estimates by dividing the sum of the product of depth and number of reads to the total number of reads. 

> SUM( DEPTH * N_READS ) / SUM(N_READS)

Then take the target coverage and divide by the mean coverage to get the fraction to subsample.

```
$ samtools view -s <fraction_to_sub-sample> -o <out.bam> <in.bam> 
```

Check the coverage again with `bedtools genomecov`

                                                                                     166,1         Bot
