#!/usr/bin/perl
use strict; use warnings;
#my $r="/home/dantakli/ref/human_g1k_v37.fasta";
my $r="/oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/human_g1k_v37.22.fasta";
my $threads = 12;
open IN, "mixed_fq.txt";
open SH, ">batch/bwa_cmd";
my $c=1;
while(<IN>){ 

	chomp;
	my $fq = $_;
	my $fq1 = $fq."1.fq";
	my $fq2 = $fq."2.fq";
	my $o = $fq;
	$o =~ s/_$/\.bwamem\.bam/;
	$o =~ s/\.$/\.bwamem\.bam/;
	if(! -e $fq1) { warn "ERROR $fq1 DOES NOT EXIST\n"; die; }
	if(! -e $fq2) { warn "ERROR $fq2 DOES NOT EXIST\n"; die; }
	my @ff = split /\//, $fq;
	my $af = pop @ff; 
	$af =~ s/snp_sim\.1000x\.//;
	$af =~ s/\.$//;
	print $af,"\n";
	my $ofh = "/home/dantakli/somatic_sim/batch/bwa/$af\.sh";
	open OUT, ">$ofh";
	print OUT "#!/bin/sh\n";
	print OUT "bwa mem -t $threads -c 500 $r $fq1 $fq2 -R \"\@RG\\tID:$af\\tSM:$af\\tLB:$af\\tPL:ILLUMINA\" | samtools view -bh -@ $threads -O BAM -o $o\n";
	my @f = split /\//, $o;
	my $iid = pop @f;
	$iid =~ s/\.bwamem\.bam//;
	my $T = "/oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/$iid\_tmp/";
	print OUT "mkdir -p $T\n";
	my $oo = $o;
	$oo =~ s/\.bam/\.sorted\.bam/;
	print OUT "samtools sort -@ $threads -o $oo -O BAM -T $T$iid $o\n";
	print OUT "rm -rf $T\n";
	print OUT "samtools index $oo\n";
	print OUT "bedtools genomecov -ibam $oo -g /home/dantakli/somatic_sim/human.hg19.genome >/home/dantakli/somatic_sim/cov/sim_snp.1000x.$af\.genomecov.txt\n";  
	close OUT;
	print SH "bash $ofh\n";
} close IN;
close SH; 
