#!/usr/bin/perl -w
use strict;

########################################################################
# gff2bed_features.pl
#
#my own gff2bed#
#June 19 2012, modified July 20 2012#
#modified 06-23-2015 to make generic

#modified 11-10-2022 to take selected features only
#note that there is limited error-checking built in, use with caution
########################################################################

#USAGE:file name on cmd line, feature name on command line 
die "USAGE: ./gff2bd gff_file_name feature_type" unless @ARGV==2;

my $filename = shift;
my $feature = shift;

open(INFILE, "<$filename") or die "bad file name: $!";

while(<INFILE>){
	chomp;
	last if $_ =~ /##FASTA/;
	next if $_ =~ /^#/;
	next unless $_ =~ /\w/;
	my @F= ();
	@F=	split(/\t/);
	next unless $F[2] eq "$feature";
	
	#THIS MIGHT NEED TO BE FIXED DEPENDING ON HOW SPECIES GFF IS FORMATTED
	$F[8] =~ /ID=(.*?);/;  #note this might not be proper format other than for FlyBase
	my $a=$1;
	
	my $astart = $F[3]-1;  #correct for zero-based coordinates going from gff to bed
	my $aend = $F[4];
		
	print  "$F[0]\t$astart\t$aend\t$a\n" ;
	
} #end while	

close(INFILE);

__END__	
	

	
	

