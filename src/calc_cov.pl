#!/usr/bin/perl
use strict; use warnings;
print "IID\tMEAN_GENOME_COVERAGE\n";
foreach my $f (glob("cov/*")){
	my $freq_sum=0;
	my $prod_sum=0;
	open IN, $f;
	my @f = split /\//, $f;
	my $id = pop @f;
	$id =~ s/_genomecov\.txt//;
	$id =~ s/\.genomecov\.txt//;
	while(<IN>){
		chomp;
		my @r = split /\t/, $_;
		#next if($r[0] ne "genome");
		next if($r[0] ne "22");
		$prod_sum += ($r[1]*$r[2]);
		$freq_sum += $r[2];
	}close IN;
	my $mean_cov = $prod_sum / $freq_sum;
	print "$id\t$mean_cov\n"; 
}
