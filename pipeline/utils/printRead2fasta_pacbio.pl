#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use IO::Compress::Gzip;

###################################################################################################################
#                                                                                                                 #
# This script serves to export PacBio reads from the original bam file and export reads equal or longer than 10kb #
# into a single fasta file.                                                                                       #
# During this process flag, chromosome and genomic location are appened to the original read name.                #
#                                                                                                                 # 
###################################################################################################################
#author: David Porubsky

my $usage = "Usage $0 <bam_1, bam_2...>\n";
die $usage if @ARGV < 1;

my @files = ();
foreach (@ARGV) { push @files, $_; }

foreach my $bam (<@files>) {
	warn("Working on $bam\n");
	my $filename = $bam;
	$filename =~ /(\w+).(\w+)\.bam/;
	$filename =  $1."_".$2;

	$filename = $filename."_PBreads.fasta";

	open (BAM, "samtools view $bam |") or die "Data processing failed";

	while(<BAM>) {
		chomp;
		my ($readName,$flag,$chrom,$pos,$mapq,$seq) = (split "\t", $_)[0,1,2,3,4,9];

		if ($chrom eq "*") {
			$chrom = 'unknown';
		}

		#filter PacBio reads
		next if $flag & 256; #skip secondary alignments	
		#next if $mapq < 254; 
		next if length($seq) <= 10000; #take only reads equal or longer than 10kb
		
		#next if $flag == 16;
		if ($flag == 16) {
			$seq =~ tr/ACGTacgt/TGCAtgca/;		
			$seq = reverse($seq);
		}
	
		print ">$readName"."_"."$flag"."_"."$chrom"."_"."$pos\n";
		print "$seq\n";
	}
}
