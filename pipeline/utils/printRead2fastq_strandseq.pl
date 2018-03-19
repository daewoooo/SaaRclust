#!/usr/bin/perl -w
use strict;
use List::Util qw/shuffle/;

##################################################################################################################
#                                                                                                                #
# This script serves to export Strand-seq reads from original bam files and export them as a single fastq file.  #
# During this process library name, flag, chromosome and genomic location are appened to the original read name. #
#                                                                                                                # 
##################################################################################################################
#author: David Porubsky

my $usage = "Usage $0 <bam_1, bam_2...>\n";
die $usage if @ARGV < 1;

my @files = ();
foreach (@ARGV) { push @files, $_; }

open OUT,'>', "strandS_libs_NA12878_list.txt" or die "Could not write to output file: $!\n";

foreach my $file (@files) {
        my $file_id = $file;

        $file_id =~ /\/(.+_\d+)_/;
        $file_id = $1;
        warn("Working on library $file_id\n");

        print OUT "$file_id\n";

		  open (BAM, "samtools view $file|") or die "Data processing failed";

		  while(<BAM>) {
			  chomp;
			  my ($readName, $flag, $chr, $pos, $seq, $qual) = (split "\t", $_)[0,1,2,3,9,10];
		
			  if ($chr eq "*") {
                $chr = 'unknown';
           }

			  #optional filtering criteria
			  #next if $flag&1024;
			  
			  if ($flag&64) {  #check if first in pair
		
			  if ($flag&16) { #read reverse strand
				  $seq =~ tr/ACGTacgt/TGCAtgca/;		
				  $seq = reverse($seq);
			  }

			  my $ID = "@".$readName."_".$file_id."_".$flag."_".$chr."_".$pos;		

			  print "$ID\n";
			  print "$seq\n";
			  print "+\n";
			  print "$qual\n";
		}
	}
}
close OUT;
