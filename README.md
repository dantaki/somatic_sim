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
9. GATK HaplotypeCaller

---


1. Generate random SNPs

```
$ randomBed -l 1 -n 1000000 -g human.hg19.22.genome  | cut -f 1-3  | sortBed >sim_snps.22.bed
```

2 Remove SNPs within excluded regions

## Excluded Regions:

* RepeatMasker 
* Short Tandem Repeats (STR)
* Segmental Duplications (SegDups)
* Unmappable Regions

```
$ intersectBed -a sim_snps.22.bed -b /home/dantakli/resources/repeatmasker/repeatMasker_hg19.bed.gz \
                                     /home/dantakli/resources/hg19_str.bed.gz \
                                     /home/dantakli/resources/hg19_segdup.bed \
                                     /home/dantakli/resources/hg19_umap_k24_mask.bed.gz \
                                     -wa -v >sim_snps.22.filtered.bed
```

3. Randomly select 10,000 SNPs

```
$ shuf sim_snps.22.filtered.bed | head -n 10000 | sortBed >sim_snps.22.filtered.10k.bed
```

4. Make a VCF with ALT alleles

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

5. GATK FastaAlternateReferenceMaker

Use GATK to create a FASTA file containing ALT alleles. This FASTA file will be used to simulate Illumina reads.

First generate a FASTA file of the contig of interest (here chromosome 22)

```
$ samtools faidx human_g1k_v37.fasta 22 >human_g1k_v37.22.fasta
$ samtools faidx human_g1k_v37.22.fasta # index the FASTA file
$ bwa index human_g1k_v37.22.fasta # bwa index the FASTA file for alignment steps
```

Generate the ALT FASTA file:

* Edit this script to define output file * 
```
bash scripts/gatk_fasta_alt_maker.sh in.vcf
```

**GATK command:**

```

java -Xmx16g -jar GenomeAnalysisTK.jar \
-T FastaAlternateReferenceMaker -R $fasta -o /oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/human_g1k_v37.22.ALT.fasta -L 22 -V $1
```


5. ART Illumina Simulator 1000x
  * REF
  * ALT

6. pysim mixture_v1.py

7. bwa mem align

8. samtools sort & sub-sample to 500x
                                                                                     166,1         Bot

