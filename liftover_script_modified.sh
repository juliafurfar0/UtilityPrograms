#!/bin/bash

####################################################################
#script for liftOver to get sequences for SCRMshaw training sets   #
# (c) Marc S. Halfon November 28, 2017                             #
#Takes a BED file of CRMs and gets aligned sequences from other    #
#Drosophila species. To do this must first convert to Dm3          #
#coordinates. This is done using liftOver. Other species           #
#coordinates are then obtained also using liftOver.                #
#The log file lists how many unmapped sequences from each.         #
#Coordinates for other species are fed to bedtools getfasta to get #
#actual FASTA files.                                               #
#Finally, the FASTA headers are modified to contain the species    #
#and genome build.                                                 #
# **NOTE that the chain file names and header modifications are    #
#hard-coded and need to be changed when changing species and/or    #
#builds.                                                           #
#                                                                  #
#revision history:                                                 #
#version 1: November 28 2017                                       #
#v2: changing name to name+ grep updated                            #

#v3: Hasiba Asma                                                   #
#changes made to run new chain files (which allow direct conversion #
#from dm6 to related species)                                       #
####################################################################


#Usage: cmd line provide (1) Dm6 bed file (2) path to liftOver (3) path to chain files (4) path to genome fastas

#check number of arguments
if [ "$#" -ne 4 ]; then
	echo "Usage: $0 Dm6_bed_file path_to_liftOver path_to_chain_files path_to_genome_fastas"
	exit 1
fi	

Dm6BED=$1
LIFTOVERDIR=$2
CHAINDIR=$3
FASTADIR=$4


#check file and directory types
if [ ! -f "$Dm6BED" ]; then
	echo "$Dm6BED not a file"
	exit 1
fi

if [ ! -d "$LIFTOVERDIR" ]; then
	echo "$LIFTOVERDIR not a valid directory"
	exit 1
fi	

if [ ! -d "$CHAINDIR" ]; then
	echo "$CHAINDIR not a valid directory"
	exit 1
fi

if [ ! -d "$FASTADIR" ]; then
	echo "$FASTADIR not a valid directory"
	exit 1
fi

# echo $LIFTOVERDIR
# echo $CHAINDIR
# echo $FASTADIR


date=`date +%m%d%Y`

#create output directory
#outdir="output.$date"
#mkdir $outdir

#create log file
logfile="log.liftover.$date"
echo "liftOver log $date" 1>$logfile

echo "Input file: $Dm6BED" 1>>$logfile

totalstarting=`wc -l $Dm6BED`
echo "number of input sequences: $totalstarting" 1>>$logfile

#----------------------------------------#
#convert Dm6 coordinates BED file to Dm3 to use with liftOver
# DM3BED="dm3.$Dm6BED"
# #echo $DM3BED
# 
# $LIFTOVERDIR/liftOver -minMatch=0.25 $Dm6BED $CHAINDIR/dm6ToDm3.over.chain $DM3BED $DM3BED.unmapped  ##1>>$logfile
# 
# NUM=`grep -vc '#' dm3.$Dm6BED.unmapped`
# echo "Dm3 unmapped: $NUM" 1>>$logfile
# 
# totalconverted=`wc -l $DM3BED`
# echo "number of dm3 sequences: $totalconverted" 1>>$logfile

#===================================================#
#species loop:
# species: Dana Dere Dgri Dmoj Dper Dpse Dsec Dsim Dvir Dyak
#this section loops through the species and performs both liftOver and getfasta
#to add/delete species, need to:
# (1) add/delete from "FOR" statement and
# (2) add/delete "CASE" with appropriate file names
# (3) add/delete header correction in following section

for SPECIES in Dana Dere Dgri Dmoj Dper Dpse Dsec Dsim Dvir Dyak
do

#----------------------------------------#
#case to set species-specific variables
case $SPECIES in
Dana)
	chainfile="dm6.GCF_000005115.1_dana_caf1.over.chain"
	fastafile="dana_caf1_genomic.fna"
	;;
Dere)
	chainfile="dm6.GCF_003286155.1_DereRS2.over.chain"
	fastafile="DereRS2_genomic.fna"
	;;
Dgri)
	chainfile="dm6.GCF_000005155.2_dgri_caf1.over.chain"
	fastafile="dgri_caf1_genomic.fna"
	;;
Dmoj)
	chainfile="dm6.GCF_000005175.2_dmoj_caf1.over.chain"
	fastafile="dmoj_caf1_genomic.fna"
	;;
Dper)
	chainfile="dm6.GCF_003286085.1_DperRS2.over.chain "
	fastafile="DperRS2_genomic.fna"
	;;
Dpse)
	chainfile="dm6.GCF_000001765.3_Dpse_3.0.over.chain"
	fastafile="Dpse_3.0_genomic.fna"
	;;
Dsec)
	chainfile="dm6.GCF_000005215.3_dsec_caf1.over.chain"
	fastafile="dsec_caf1_genomic.fna"
	;;
Dsim)
	chainfile="dm6.GCF_000754195.2.dsim.over.chain"
	fastafile="dsim_genomic.fna"
	;;
Dvir)
	chainfile="dm6.GCF_000005245.1_dvir_caf1.over.chain"
	fastafile="dvir_caf1_genomic.fna"
	;;
Dyak)
	chainfile="dm6.GCF_000005975.2_dyak_caf1.over.chain"
	fastafile="dyak_caf1_genomic.fna"
	;;
esac

#----------------------------------------#								
#liftover and bedtools 'getfasta'
echo Starting $SPECIES

#$LIFTOVERDIR/liftOver -minMatch=0.25 $DM3BED $CHAINDIR/$chainfile $SPECIES.$DM3BED $SPECIES.$DM3BED.unmapped 
# echo $Dm6BED 
# echo $CHAINDIR/$chainfile
# echo $SPECIES.$Dm6BED
# echo $SPECIES.$Dm6BED.unmapped
# echo Starting
echo $LIFTOVERDIR/liftOver -minMatch=0.25 $Dm6BED $CHAINDIR/$chainfile $SPECIES.$Dm6BED $SPECIES.$Dm6BED.unmapped

$LIFTOVERDIR/liftOver -minMatch=0.25 $Dm6BED $CHAINDIR/$chainfile $SPECIES.$Dm6BED $SPECIES.$Dm6BED.unmapped 
echo ""

#NUM=`grep -vc '#' $SPECIES.$DM3BED.unmapped`
NUM=`grep -vc '#' $SPECIES.$Dm6BED.unmapped`
echo "$SPECIES unmapped: $NUM" 1>>$logfile

#bedtools getfasta -name+ -fi $FASTADIR/$fastafile -bed $SPECIES.$DM3BED -fo $SPECIES.$Dm6BED.fa
bedtools getfasta -name+ -fi $FASTADIR/$fastafile -bed $SPECIES.$Dm6BED -fo $SPECIES.$Dm6BED.fa

done   #end of loop

##
bedtools getfasta -name+ -fi $FASTADIR/dm6.fa -bed $Dm6BED -fo Dmel.$Dm6BED.fa


##
#===================================================#
#header correction
#this section puts the species name and build into the FASTA header
#must check proper syntax for each species/build

# perl -i -pe 's/\:scaffold/Dana3\:scaffold/' Dana.$Dm6BED.fa
# perl -i -pe 's/\:scaffold/Dere2\:scaffold/' Dere.$Dm6BED.fa
# perl -i -pe 's/\:scaffold/Dgri2\:scaffold/' Dgri.$Dm6BED.fa
# perl -i -pe 's/\:scaffold/Dmoj3\:scaffold/' Dmoj.$Dm6BED.fa
# perl -i -pe 's/\:super/Dper\:super/' Dper.$Dm6BED.fa
# perl -i -pe 's/\:chr/Dpse4\:chr/' Dpse.$Dm6BED.fa
# perl -i -pe 's/\:super/DSec1\:super/' Dsec.$Dm6BED.fa
# perl -i -pe 's/\:chr/DSim1\:chr/' Dsim.$Dm6BED.fa
# perl -i -pe 's/\:scaffold/Dvir3\:scaffold/' Dvir.$Dm6BED.fa
# perl -i -pe 's/\:chr/Dyak2\:chr/' Dyak.$Dm6BED.fa

# 
perl -i -pe 's/\:scaffold/Dana3\:scaffold/' Dana.$Dm6BED.fa
perl -i -pe 's/\:NW/Dere2\:NW/' Dere.$Dm6BED.fa
perl -i -pe 's/\:NW/Dgri2\:NW/' Dgri.$Dm6BED.fa
perl -i -pe 's/\:NW/Dmoj3\:NW/' Dmoj.$Dm6BED.fa
perl -i -pe 's/\:NW/Dper\:NW/' Dper.$Dm6BED.fa
perl -i -pe 's/\:NW/Dpse4\:NW/' Dpse.$Dm6BED.fa
perl -i -pe 's/\:NC/Dpse4\:NC/' Dpse.$Dm6BED.fa
perl -i -pe 's/\:NW/DSec1\:NW/' Dsec.$Dm6BED.fa
perl -i -pe 's/\:NT/DSim1\:NT/' Dsim.$Dm6BED.fa
perl -i -pe 's/\:NC/DSim1\:NC/' Dsim.$Dm6BED.fa
perl -i -pe 's/\:NW/Dvir3\:NW/' Dvir.$Dm6BED.fa
perl -i -pe 's/\:NT/Dyak2\:NT/' Dyak.$Dm6BED.fa
perl -i -pe 's/\:NC/Dyak2\:NC/' Dyak.$Dm6BED.fa
perl -i -pe 's/\:chr/Dmel\:chr/' Dmel.$Dm6BED.fa


################################
#notes:
# can use grep to extract desired subset of FASTA sequences. Code is:
# grep -h -A 1 '[name(s) of CRM]' *.fa > outfile
# perl -i -pe 's/^--\n//' outfile




