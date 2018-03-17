#!/usr/bin/perl -w
use strict;
use Data::Dumper;

#my $inputFile = shift;
#open IN,'<', $inputFile or die "Could not read the file: $!\n";

print "SSreadNames\tSSlibNames\tSSflag\tSSchrom\tSSpos\tstrand\tPBreadNames\tPBflag\tPBchrom\tPBpos\tPBreadLen\tTargetCoordStart\tTargetCoordend\tMatchedBasesWithGaps\n";

while (<STDIN>) {
	chomp;
	my($SSreadNames, $strand, $PBreadNames, $PBreadLen, $TargetCoordStart, $TargetCoordend, $MatchedBasesWithGaps) = (split "\t", $_)[0,4,5,6,7,8,10];

	#process StrandS read names
	my @SSreadNames_fields = (split "_", $SSreadNames);
	$SSreadNames = $SSreadNames_fields[0];	
	my $SSlibNames = join("_", ($SSreadNames_fields[1],$SSreadNames_fields[2]));
	my $SSflag = $SSreadNames_fields[3];
	my $SSchrom = $SSreadNames_fields[4];
	my $SSpos = $SSreadNames_fields[5];

	next if $SSflag&1024;

	#process PacBio read names
	my @PBreadNames_fields = (split "_", $PBreadNames);
	$PBreadNames = join("_", @PBreadNames_fields[0..5]);
	my $PBflag = $PBreadNames_fields[7];
	my $PBchrom = $PBreadNames_fields[8];
	my $PBpos = $PBreadNames_fields[9];

        print "$SSreadNames\t$SSlibNames\t$SSflag\t$SSchrom\t$SSpos\t$strand\t$PBreadNames\t$PBflag\t$PBchrom\t$PBpos\t$PBreadLen\t$TargetCoordStart\t$TargetCoordend\t$MatchedBasesWithGaps\n";
}
