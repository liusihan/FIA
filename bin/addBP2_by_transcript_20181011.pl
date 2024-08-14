#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
#use File::Basename;
use Getopt::Long;

#############################
#   2018-10-11              #
#	Yuanyuan Zeng  			#
#############################


my ($annofile,$transcriptfile, $linkagefile,$outputfile,$check);

GetOptions(
	'anno=s' => \$annofile,
	't=s' => \$transcriptfile,
	'l=s' => \$linkagefile,
	'o=s' => \$outputfile,
	'check=s' => \$check,
);

die "Usage: perl $0 -anno vep.vcf -t transcript_file -l linkage_file -o output.vcf \n" unless (defined($annofile) && defined($outputfile));

my %refseq;


open(REFSEQ, $transcriptfile);
open(LINK, $linkagefile);
if($annofile =~ /\.gz$/){
	open(VEP, "gzip -dc $annofile|");
}
else{
	open(VEP, $annofile);
}

open(OUT, ">", $outputfile);
open(LOG, ">", "$outputfile.log");



my %hash_P;
#my %hash_V;
my %all_var;
my %stat_pathogenic;

while(<REFSEQ>){
	chomp;
	if(/chrom/){next;}
	my @line = split;
	my $id = join("_", @line[0..3]);
#	$all_var{$id} = 1;
	#$stat_pathogenic{$line[4]}++;
#	if($line[4] eq "P" || $line[4] eq "PW") {$hash_P{$id} = 1;}
	$hash_P{$id} = 1;
#	if($line[4] eq "VW") {$hash_V{$id} = 1;}
}
close REFSEQ;
#print LOG Dumper(\%hash_P),"\n";

my %link_variants;
my %hash_bp2;
while (<LINK>){
	chomp;
	if( /^sampleid/ ){ next; }
	my @line = split(/\t/,$_);
	my $id = join("_", @line[1..4]);
	push @{$link_variants{$line[0]}{$line[1]}}, $id;
}
close LINK;
#print LOG Dumper(\%link_variants),"\n";


foreach my $s(keys %link_variants){
	foreach my $chr(keys %{$link_variants{$s}}){
		my $flag = 0;
		my @array_V;
		foreach my $v(@{$link_variants{$s}{$chr}}){
			if(exists($hash_P{$v}))
			{
				$flag++;
			}
			else{
				push @array_V, $v;
			}
		}
		if( $flag > 0 && @array_V > 0 ){
			$hash_bp2{$_} = 1 foreach @array_V;
		}
	}
}
#close LOG;
print LOG Dumper(\%hash_bp2),"\n";
#=====================================================

my @title;
my @t_index;
#my %vep;
my %stat;
$stat{"*"} = 0;
$stat{"BP2"} = 0;
$stat{"var"} = 0;
my %p_stat;


## read vep result
while(<VEP>){
	chomp;
	if(/^##/){next;}
	if(/^#CHROM/){
		@title[0..8] = split(/\t/,$_);
		print OUT join("\t","chrom","pos","ref","alt","transcript", "BP2"),"\n";
		next;
	}
	$stat{"total"}++;
	my @line;
	@line[0..8] = split(/\t/, $_);
	if($line[4] eq "*"){
		$stat{"*"}++;
		next;
	}
#===========================
	my $id = join("_", @line[0..1], @line[3..4]);
	my $bp2_note = "NA";
	my $p_note = "Pathogenic_NotRecord";
	if( exists($all_var{$id}) )
	{ 
		$p_note = $all_var{$id}; 
	}

	if( exists($hash_bp2{$id}) ) 
	{ 
		$bp2_note = 3;
		$stat{"BP2"}++;
		$p_stat{$p_note}++;	
		$stat{"var"}++;
	}

	my $feature;
	if($line[7] =~ /(.*);VEP_Feature=(\S+);VEP_BIOTYPE=(.*)/){$feature = $2 . "|end";}
	if(!$feature){next;}

	my @feature_all = split(/\|/,$feature);
	foreach my $i(0..($#feature_all - 1)){

#		$feature_all[$i] =~ s/\.\d+//g;
		print OUT join("\t", @line[0..1], @line[3..4], $feature_all[$i], $bp2_note),"\n";
	}
}

my @p = keys %p_stat;

print LOG join("\t", "total", "ALT=*","BP2_variants","BP2_transcripts", @p),"\n";
print LOG join("\t", $stat{"total"}, $stat{"*"}, $stat{"var"}, $stat{"BP2"});
print LOG "\t$p_stat{$_}" foreach @p;
print LOG "\n";

close VEP;
close OUT;
close LOG;


