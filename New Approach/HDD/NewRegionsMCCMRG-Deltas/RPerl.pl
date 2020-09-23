#!/usr/bin/perl
use warnings;
use strict;

my @MCC = (1..40);
my $MCCInterv = 0.5;
my @MRG = (1..10);
my $MRGInterv = 10;
my @MCCRange = map {$_ * $MCCInterv} @MCC;
my @MRGRange = map {$_ * $MRGInterv} @MRG;
my $MCCRangeLength = @MCCRange;
my $MRGRangeLength = @MRGRange;
my $MCC_MRG = $MCCRangeLength * $MRGRangeLength;
my @MCC_MRG = (1..$MCC_MRG);
# my $output = `sudo Rscript ./NewRegionsMCCMRG-Deltas.R ; free -h && sudo sysctl -w vm.drop_caches=3 && sudo sync && echo 3 | sudo tee /proc/sys/vm/drop_caches && free -h`;

foreach(@MCC_MRG){
	system('sudo Rscript ./NewRegionsMCCMRG-Deltas.R ; free -h && sudo sysctl -w vm.drop_caches=3 && sudo sync && echo 3 | sudo tee /proc/sys/vm/drop_caches && free -h');
	# $output;
}

#print("$MCC\n");
#print 10 * 10
#MRG <- c(seq(10,100, 10))

#sudo Rscript ./RequiredPackages.R 
# my @r = (1..390);

# for(@r){
# 	print("$_","\n");
# }

# for (my $r = 1; $r <= 390;$r += 5) {print "$range\n"; }
