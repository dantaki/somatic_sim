#!/usr/bin/perl
use strict; use warnings;
# output directory
my $dir = "/oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/split_mix/";

# reference FASTQ prefix
my $ref = "/oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/split/ref.22.1000x."; 

# alt FASTQ prefix
my $alt = "/oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/split/alt.22.1000x.";

# split FASTQ suffixes (if you split the FASTQ)
my @a = qw/ 00. 01. 02. 03. 04./;

# file of commands to submit as a job-array
open SH, ">batch/mix_cmd";

# allele fractions to test
my @params = qw/0.50 0.25 0.20 0.15 0.10 0.05 0.04 0.03 0.02 0.01/;

# foreach split FASTQ
foreach my $f (@a) { 
	# $altm is the mixture coefficient for the ALT FASTQ
	foreach my $altm (@params){
		my $refm = 1-$altm;
		# output prefix for the .conf file
		my $conf = "/home/dantakli/somatic_sim/snp_mix_confs/$f\.$altm\.conf";
		open OUT, ">$conf";
		# define FASTQ prefixes
		my $r = $ref.$f;
		my $a = $alt.$f;
		
		print OUT "$r\t$refm\n";
		print OUT "$a\t$altm\n"; 
		close OUT;

		# output FASTQ filename
		my $o_prefix= "sim_snps";
		my $o = $f; 
		$o =~ s/\.//;
		my $oaltm= $altm;
		$oaltm =~ s/^0\.//;
		
		print SH "python /home/dantakli/somatic_sim/pysim/mixture_v1.py -i $conf -o $dir$o_prefix\.$oaltm"."AF.$o\n"; 
	}
}
close SH;
