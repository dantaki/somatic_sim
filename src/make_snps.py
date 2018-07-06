#!/usr/bin/env python2
import pysam
from random import shuffle
import sys
fasta = pysam.FastaFile("/home/dantakli/ref/human_g1k_v37.fasta")
alts = ['A','T','G','C']
ofh = open('sim_snps.22.filtered.vcf','w')
ofh.write('##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
c=1
with open(sys.argv[1]) as f:
	for l in f:
		vid = 'SIMSNP{}'.format(c) 
		c+=1
		chrom, start, end = l.rstrip().split('\t')
		ref = fasta.fetch(region='{}:{}-{}'.format(chrom,end,end))
		tmp_alt = [x for x in alts if x != ref]
		shuffle(tmp_alt)
		ofh.write('{}\t{}\t{}\t{}\t{}\t.\t.\t.\n'.format(chrom,end,vid,ref,tmp_alt[0]))
ofh.close()
