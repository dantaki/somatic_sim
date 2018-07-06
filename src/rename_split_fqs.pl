#!/usr/bin/perl
use strict; use warnings;
#foreach my $f(glob("/oasis/projects/nsf/ddp195/dantakli/somatic_sim/*_split/*")){
foreach my $f(glob("/oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/split/*")){
	next if($f =~ /fq$/);
	my @f = split /\//, $f;
	my $fq = pop @f;
	my $dir = join "/", @f;
	my @i = split /\./, $fq;
	my $d = pop @i;
	my $fqnm = pop @i;
	my $final = join ".", @i;
	$final = $dir."/".$final.".$d.$fqnm.fq";
	print "mv $f $final\n"; 	
}
