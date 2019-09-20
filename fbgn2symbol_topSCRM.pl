#!/usr/bin/perl -w
use strict;

##############################
# script to put gene symbols in one of the FBgn fields in results
# file from "Generate_topN_SCRMhits" v3.4 output
#
# (c) MSH  October 26, 2017
# revised 09-2019 to allow multiple FBgn's per field
##############################
# requires input files:
# (1) topN hits file
# (2) flybase conversion file 'flybase_nameconversion'
# latter download from FlyBase FBgn <=> Annotation ID

#use only primary IDs and symbols, too confusing otherwise

my %fbgns;

my ($hitsfile, $IDfile) = @ARGV;

#load and hash the flybase data
open(ID, "<$IDfile") or die "Can't open ID file $!";
while(<ID>){
	chomp;
	my @line = split;
	$fbgns{$line[2]} = $line[0];
}
close(ID);

open(HITS, "<$hitsfile") or die "Can't open hits file $!";	
while(<HITS>){
	chomp;
	my @line = split;
	
		my @sym = ();
		foreach my $id (split(/,/, $line[5])){
			
	
			if(exists $fbgns{$id}){
				push(@sym, $fbgns{$id});
			}	else {
				push(@sym, $id);
			}	
		}
		$line[5] = join(',', @sym);
	

	
		@sym = ();
			foreach my $id (split(/,/, $line[10])){
			
				
				if(exists $fbgns{$id}){
					push(@sym, $fbgns{$id});
				}	else {
					push(@sym, $id);
				}	
			}
			$line[10] = join(',', @sym);
		
	print join("\t", @line);
	#print "$line[5] and $line[10]";
	print "\n";
}
close(HITS);	