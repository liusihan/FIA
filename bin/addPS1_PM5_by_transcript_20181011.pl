#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
#use File::Basename;
use Getopt::Long;

#############################
#   2018-10-12              #
#	Yuanyuan Zeng  			#
#############################


my ($annofile,$transcriptfile,$outputfile,$check);

GetOptions(
	'anno=s' => \$annofile,
	't=s' => \$transcriptfile,
	'o=s' => \$outputfile,
	'check=s' => \$check,
);

die "Usage: perl $0 -anno vep.vcf -t transcript_file -o output.vcf \n" unless (defined($annofile) && defined($outputfile));

my %refseq;

open(REFSEQ, $transcriptfile);
if($annofile =~ /\.gz$/){
	open(VEP, "gzip -dc $annofile|");
}
else{
	open(VEP, $annofile);
}
open(OUT, ">", $outputfile);
open(LOG, ">", "$outputfile.log");
my %ps1_id;
my %pm5_id;
my %p_variants;
my @title_indx;
##===============  part 1 ===================================================
while(<REFSEQ>){
	chomp;
	if(/VEP_HGVSp/)
	{
		my @title = split(/\t/, $_);
		$title_indx[0] = (grep{$title[$_] eq "VEP_HGVSp"} 0..$#title)[0];
		$title_indx[1] = (grep{$title[$_] eq "VEP_Consequence"} 0..$#title)[0];
		#print $title_indx[1],"\n";
		next;
	}
	my @line = split(/\t/,$_);
	if(@line < 7){next;}
	my $varid = join("_", @line[0..3]);
	$p_variants{$varid} = 1;
	my $hgvsp = $line[$title_indx[0]];
	my $consequences = $line[$title_indx[1]];
#	print Dumper(\@line),"\n";
#	print $consequences,"\n";

#	my @proteins;
#	if(defined $hgvsp){
	my @proteins = split(/\|/, $hgvsp);
	# }
	# else
	# {
	# next;
	# }
	my @consequence = split(/\|/, $consequences);

	foreach my $i(0..$#proteins){
		my $p = $proteins[$i];
		if($p eq ""){next;}
		#if($consequence[$i] =~ /missense_variant/){
			$ps1_id{$p}++;
			if($consequence[$i] =~ /missense_variant/){
				$p =~ s/[a-zA-Z]+$//;
				$pm5_id{$p}++;
			}
	}
}
close REFSEQ;

#print Dumper(\%pm5_id),"\n";

#===============  part 2 ============================================

my @title;
my @t_index;
my %stat;
$stat{"*"} = 0;
$stat{"PS1"} = 0;
$stat{"PM5"} = 0;
## read vep result
while(<VEP>){
	chomp;
	if(/^##/){next;}
	if(/^#CHROM/){
		print OUT join("\t", "chrom","pos","ref","alt","transcript", "PS1","PM5"), "\n";
		next;
	}
	$stat{"total"}++;
	my @line;
	@line[0..8] = split(/\t/, $_);
	if($line[4] eq "*"){
		$stat{"*"}++;
		next;
	}

	my $hgvsp;
	my $consequence;
	my $feature;

	if($line[7] =~ /(.*);VEP_HGVSp=(\S+);VEP_cDNA_position=(.*)/){$hgvsp = $2 . "|end";}
	if($line[7] =~ /(.*);VEP_Consequence=(\S+);VEP_IMPACT=(.*)/){$consequence = $2 . "|end";}
	if($line[7] =~ /(.*);VEP_Feature=(\S+);VEP_BIOTYPE=(.*)/){$feature = $2 . "|end";}
	if(!$feature){next;}

	my @feature_all = split(/\|/,$feature);

	if(!defined($hgvsp))
	{
#		print "$feature\n";
#		print Dumper(\@feature_all), "\n";

		foreach my $i(0..($#feature_all - 1)){
#			$feature_all[$i] =~ s/\.\d+//g;
#			$i =~ s/\.\d+//g;
			print OUT join("\t", @line[0..1], @line[3..4], $feature_all[$i], "NA", "NA"), "\n";
		}
		next;
	}
#	print Dumper(\@feature_all), "\n";


	my @hgvsp_all = split(/\|/,$hgvsp);
	my @consequence_all = split(/\|/,$consequence);
	my $id = join("_", @line[0..1], @line[3..4]);

	my $varid = join("_", @line[0..1], @line[3..4]);
	foreach my $i(0..($#feature_all - 1)){

#		$feature_all[$i] =~ s/\.\d+//g;

		if(exists($p_variants{$varid}))
		{
			print OUT join("\t", @line[0..1], @line[3..4], $feature_all[$i], "NA","NA"), "\n";
			next;
		}

		my $p = $hgvsp_all[$i];
		$p =~ s/%3D/=/;
		if (exists($ps1_id{$p}))
		{
			print OUT join("\t", @line[0..1], @line[3..4], $feature_all[$i], "1","NA"), "\n";
			$stat{"PS1"}++;
			$stat{"var_ps1"}{$id}++;
			next;
		}
		if($consequence_all[$i] =~ /missense_variant/)
		{
			$p =~ s/[a-zA-Z]+$//;
			if(exists($pm5_id{$p}))
			{ 
				if($pm5_id{$p}==1)
					{print OUT join("\t", @line[0..1], @line[3..4], $feature_all[$i], "NA","2"), "\n";}
				if($pm5_id{$p}>=2)
					{print OUT join("\t", @line[0..1], @line[3..4], $feature_all[$i], "NA","1"), "\n";}
				$stat{"PM5"}++;
				$stat{"var_pm5"}{$id}++;
				next;
			}
		}
		print OUT join("\t", @line[0..1], @line[3..4], $feature_all[$i], "NA","NA"), "\n";
		
	}
}

my @n_ps1 = keys %{$stat{var_ps1}};
my @n_pm5 = keys %{$stat{var_pm5}};
print LOG join("\t", "total", "ALT=*","PS1_Variants", "PM5_Variants", "PS1_Transcripts", "PM5_Transcripts"),"\n";
print LOG join("\t", $stat{"total"}, $stat{"*"}, @n_ps1 + 0, @n_pm5 + 0, $stat{"PS1"}, $stat{"PM5"}), "\n";
close VEP;
close OUT;
close LOG;
