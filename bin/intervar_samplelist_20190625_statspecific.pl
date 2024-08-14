#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use List::Util qw/sum/;
use File::Basename;
use Storable qw(dclone);
use Cwd 'abs_path';
use File::Spec;

die "perl $0 <ar_file> <ad_file> <sample_list> <manual_evidence_file> <outdir>\n" unless(@ARGV == 5);
my $ar_file = $ARGV[0];
my $ad_file = $ARGV[1];
my $samplelist = $ARGV[2];
my $manual_evidence = $ARGV[3];
my $outdir = $ARGV[4];

# opendir(DIR, $sample_dir);
# my @samplelist = readdir DIR;
## sample list
open(FILE, $samplelist) or die $!;
my $content = do {local $/=undef;<FILE>;};
my @samplelist = split(/\n/, $content);
#print Dumper(\@samplelist);
##
my %db_ar;
my %db_ad;
my %pb;

my $absolute_dir = abs_path($manual_evidence);
my $evidence_mode = File::Spec->catfile(dirname($absolute_dir), "evidence_genetic_mode.txt");
my $B_P_file = File::Spec->catfile(dirname($absolute_dir),"B_P_conflict_sites.txt");

open(my $IN, '<', $evidence_mode);
while(<$IN>){
	chomp;
	my @line = split(/\t/, $_);
	$pb{$line[1]} = $line[3];
}
close ($IN);

my %special_list;
open(my $IN1, '<', $B_P_file);
while (<$IN1>) {
	chomp;
	my @line = split(/\t/,$_);
	my $varid = join("_", @line[0..3]);
	my $gene = $line[4];
	$special_list{$gene}{$varid} = 1;
}
close ($IN1);

#print Dumper(\%pb), "\n";

open(AR, $ar_file);
open(AD, $ad_file);

my @title_ar;
my @title_ad;

while(my $line1 = <AR>){
	my $line2 = <AD>;
	chomp $line1;
	chomp $line2;
	if($line1 =~ /^chrom/){
		@title_ar = split(/\t/, $line1);
		@title_ad = split(/\t/, $line2);
		next;
	}

	my @line_ar = split(/\t/, $line1);
	my @line_ad = split(/\t/, $line2);
	my $id = join("_", @line_ar[0..4]);

## BP2,BS2 conflict ##############
my $bp2_ar = (grep{$title_ad[$_] eq "BP2"} 5..$#title_ad)[0];
my $bs2_ar = (grep{$title_ad[$_] eq "BS2"} 5..$#title_ad)[0];
if(($line_ad[$bs2_ar] ne "NA") && ($line_ad[$bp2_ar] ne "NA")){
	$line_ad[$bs2_ar] = "NA";
}

## special variants ###########################
	my $id1 = join("_", @line_ar[0..3]);
	my $pvs1_ar = (grep{$title_ar[$_] eq "PVS1"} 5..$#title_ar)[0];
	my $pvs1_ad = (grep{$title_ad[$_] eq "PVS1"} 5..$#title_ad)[0];
	my $pp3_ar = (grep{$title_ar[$_] eq "PP3"} 5..$#title_ar)[0];
	my $pp3_ad = (grep{$title_ad[$_] eq "PP3"} 5..$#title_ad)[0];
	if($id1 eq "7_107302082_A_T" || $id1 eq "22_38380345_C_T" || $id1 eq "13_20766921_C_T"){
		$line_ar[$pvs1_ar] = 0;
		$line_ad[$pvs1_ad] = 0;
	}
	if($id1 eq "7_107338557_G_GT" || $id1 eq "9_75387472_G_GT" || $id1 eq "X_82764417_G_T"){
		$line_ar[$pvs1_ar] = "NA";
		$line_ad[$pvs1_ad] = "NA";
	}
	if($line_ar[$pvs1_ar] ne "NA"){$line_ar[$pp3_ar] = "NA"}
	if($line_ad[$pvs1_ad] ne "NA"){$line_ad[$pp3_ad] = "NA"}
#################################
	foreach my $i(0..3){
		${$db_ar{$id}{"P"}}[$i] = 0;
		${$db_ar{$id}{"B"}}[$i] = 0;
		${$db_ad{$id}{"P"}}[$i] = 0;
		${$db_ad{$id}{"B"}}[$i] = 0;
	}
	foreach my $n(5..$#line_ar){
		if($line_ar[$n] ne "NA"){${$db_ar{$id}{$pb{$title_ar[$n]}}}[$line_ar[$n]]++;}
	}
	foreach my $m(5..$#line_ad){
		if($line_ad[$m] ne "NA"){${$db_ad{$id}{$pb{$title_ad[$m]}}}[$line_ad[$m]]++;}
	}
	my @idx1 = grep{$line_ar[$_] ne "NA"} 5..$#line_ar;
	my @idx2 = grep{$line_ad[$_] ne "NA"} 5..$#line_ad;

	#@{$db_ar{$id}{"info"}} = ();
	#@{$db_ad{$id}{"info"}} = ();

	if(@idx1 > 0){
		foreach my $index1(@idx1){
			#my $evid = join("=", $title_ar[$index1], $line_ar[$index1]);
			#push @{$db_ar{$id}{"info"}}, $evid;
			my $evi = $title_ar[$index1];
			$db_ar{$id}{"info"}{$evi} = $line_ar[$index1];
		}
		# $db_ar{$id}{"info"} = join("|", @title_ar[@idx1]);
		# $db_ar{$id}{"value"} = join("|", @line_ar[@idx1]);
	}
	# else{
	# 	$db_ar{$id}{"info"} = "";
	# 	$db_ar{$id}{"value"} = "";
	# }
#	$db_ad{$id}{"info"} = join("=",join("|", @title_ad[@idx2]), join("|", @line_ad[@idx2]));
	if(@idx2 > 0){
		foreach my $index2(@idx2){
			#my $evid = join("=", $title_ad[$index2], $line_ad[$index2]);
			#push @{$db_ad{$id}{"info"}}, $evid;
			my $evi = $title_ad[$index2];
			$db_ad{$id}{"info"}{$evi} = $line_ad[$index2];
		}

		# $db_ad{$id}{"info"} = join("|", @title_ad[@idx2]);
		# $db_ad{$id}{"value"} = join("|", @line_ad[@idx2]);
	}
	# else{
	# 	$db_ad{$id}{"info"} = "";
	# 	$db_ad{$id}{"value"} = "";
	# }
}
close AR;
close AD;
#print Dumper(\%{$db_ar{"17_18022710_C_G"}});
### manual evidence
my %pp1;
my %ps2;

open(IN, $manual_evidence);
while(<IN>){
	chomp;
	my @line = split(/\t/, $_);
	if($line[0] eq "sampleid"){next;}
	my $id = join("_", @line[1..4]);
	my $sample = $line[0];
	if($line[5] ne "NA"){
		$pp1{$sample}{$id} = $line[5];
	}
	if($line[6] ne "NA"){
		$ps2{$sample}{$id} = $line[6];
	}
}
close IN;

#print Dumper(\%db_ar);
######################
my @all_evidence = keys %pb;
my %evi_n;
foreach my $n(0..$#all_evidence){
	$evi_n{$all_evidence[$n]} = $n;
}


my $sample_prefix;
foreach my $sample(@samplelist){
#	if($sample eq "." || $sample eq ".."){next;}
	#print $sample, "\n";
	$sample_prefix = basename($sample, ".tsv");
#	print $sample_prefix,"\n";
	open(my $H, $sample);
	open(my $OUT, ">$outdir/${sample_prefix}.classfy.tsv");
	open(my $OUT1, ">$outdir/${sample_prefix}.for_stat.AR.txt");
	open(my $OUT2, ">$outdir/${sample_prefix}.for_stat.AD.txt");

	process($H, $OUT, $OUT1, $OUT2);

	close $H;
	close $OUT;
	close $OUT1;
	close $OUT2;
}

sub process {
	my ($H, $OUT, $OUT1, $OUT2) = @_;
	my @title_sample;
	while(<$H>){
		chomp;
		my @line = split(/\t/, $_);
		my (@output_ar, @output_ad);
		foreach(0..$#all_evidence){$output_ar[$_] = "NA";}
		foreach(0..$#all_evidence){$output_ad[$_] = "NA";}
		#my @output_ar = (split //, "NA" x @title_sample);
		#my @output_ad = (split //, "NA" x @title_sample);
		#@output_ar[0..$#all_evidence] = "NA";
		#@output_ad[0..$#all_evidence] = "NA";
		my $pp4idx = $#line;
		if($line[0] eq "chrom"){
			@title_sample = @line;
			print $OUT join("\t", @title_sample[0..($#title_sample - 1)], "AR_classify","AR_P_classify", "AR_B_classify", "AR_evidence", "AD_classify","AD_P_classify", "AD_B_classify", "AD_evidence"),"\n";
			print $OUT1 join("\t", @title_sample[0..($#title_sample - 1)], "AR_classify","AR_P_classify", "AR_B_classify", @all_evidence),"\n";
			print $OUT2 join("\t", @title_sample[0..($#title_sample - 1)], "AD_classify","AD_P_classify", "AD_B_classify", @all_evidence), "\n";

			next;
		}
		my $id = join("_", @line[0..4]);
		if(!exists($db_ar{$id})){
			next;
		}

#		print $id,"\n";
		my %evidence_ar = %{dclone(\%{$db_ar{$id}})};
		my %evidence_ad = %{dclone(\%{$db_ad{$id}})};
#		if($id eq "17_18022710_C_G_ENST00000205890"){print Dumper(\%evidence_ar);}
#		if($id eq "17_18022710_C_G_ENST00000205890"){print Dumper(\%{$db_ad{$id}});}

#		if($id eq "10_55663028_A_AT_ENST00000320301") {print Dumper(\%evidence_ad);}

#		my $info_ad = join("=",$db_ad{$id}{"info"}, $db_ad{$id}{"value"});
#		my $info_ar = join("=",$db_ar{$id}{"info"}, $db_ar{$id}{"value"});
		#my $info_ad = join(";", @{$db_ad{$id}{"info"}});
		#my $info_ar = join(";", @{$db_ar{$id}{"info"}});

		if($line[$pp4idx] ne "NA"){
			${$evidence_ad{"P"}}[$line[$pp4idx]]++;
			${$evidence_ar{"P"}}[$line[$pp4idx]]++;

			$evidence_ad{"info"}{"PP4"} = $line[$pp4idx];
			$evidence_ar{"info"}{"PP4"} = $line[$pp4idx];

# 			if(@{$db_ad{$id}{"info"}} > 0){
# 				$info_ad = join(";", @{$db_ad{$id}{"info"}}, join("=", "PP4", $line[$pp4idx]));				
# #				$info_ad = join("=", join("|",$db_ad{$id}{"info"}, "PP4"), join("|", $db_ad{$id}{"value"}, $line[$pp4idx]));
# 			}
# 			else{
# 				$info_ad = "PP4=" . $line[$pp4idx];
# 			}
# 			if(@{$db_ar{$id}{"info"}} > 0 ){
# 				$info_ar = join(";", @{$db_ar{$id}{"info"}}, join("=", "PP4", $line[$pp4idx]));
# #				$info_ar = join("=", join("|",$db_ar{$id}{"info"}, "PP4"), join("|", $db_ar{$id}{"value"}, $line[$pp4idx]));
# 			}
# 			else{
# 				$info_ar = "PP4=" . $line[$pp4idx];
# 			}
		}
		# else{
		# 	if(!exists($evidence_ar{"P"})){
		# 		$evidence_ar{"P"} =[0, 0, 0, 0]
		# 		$evidence_ar{"B"} =[0, 0, 0, 0]
		# 		$evidence_ad{"P"} =[0, 0, 0, 0]
		# 		$evidence_ad{"B"} =[0, 0, 0, 0]
		# 	}
		# }
###mannual evidence
#		${$db_ar{$id}{$pb{$title_ar[$n]}}}[$line_ar[$n]]++;
		my $id_manual = join("_", @line[0..3]);
		if(exists($pp1{$sample_prefix}) && exists($pp1{$sample_prefix}{$id_manual})){
			my $strenth1 = $pp1{$sample_prefix}{$id_manual};
			$evidence_ar{"P"}[$strenth1]++;
			$evidence_ad{"P"}[$strenth1]++;
			$evidence_ad{"info"}{"PP1"} = $strenth1;
			$evidence_ar{"info"}{"PP1"} = $strenth1;

			#$info_ar = $info_ar . ";PP1=". $strenth1;
			#$info_ad = $info_ad . ";PP1=". $strenth1;
		}
		if(exists($ps2{$sample_prefix}) && exists($ps2{$sample_prefix}{$id_manual})){
			my $strenth2 = $ps2{$sample_prefix}{$id_manual};
			$evidence_ar{"P"}[$strenth2]++;
			$evidence_ad{"P"}[$strenth2]++;
			$evidence_ad{"info"}{"PS2"} = $strenth2;
			$evidence_ar{"info"}{"PS2"} = $strenth2;

			#$info_ar = $info_ar . ";PS2=". $strenth2;
			#$info_ad = $info_ad . ";PS2=". $strenth2;
		}

		my @ar_result = classfy(%evidence_ar);
		my @ad_result = classfy(%evidence_ad);
		my $info_ar = "NA";
		my $info_ad = "NA";

		if(exists($evidence_ar{"info"})){
			my %hash_ar=%{$evidence_ar{"info"}};
			
			#foreach my $e(keys %{$evidence_ar{"info"}}){
			foreach my $e(sort keys %hash_ar){
				$output_ar[$evi_n{$e}] = $evidence_ar{"info"}{$e};
				if($info_ar eq "NA"){
					$info_ar = $e . "=" . $evidence_ar{"info"}{$e};
				}
				else{
					$info_ar = $info_ar . ";" . $e . "=" . $evidence_ar{"info"}{$e};
				}
			}
		}
		#foreach my $e(keys %{$evidence_ad{"info"}}){
		if(exists($evidence_ad{"info"})){
			my %hash_ad=%{$evidence_ad{"info"}};
			foreach my $e(sort keys %hash_ad){
				$output_ad[$evi_n{$e}] = $evidence_ad{"info"}{$e};
				if($info_ad eq "NA"){
					$info_ad = $e . "=" . $evidence_ad{"info"}{$e};
				}
				else{
					$info_ad = $info_ad . ";" . $e . "=" . $evidence_ad{"info"}{$e};
				}
			}
		}
		## special list
		my $geneid = $line[5];
		if($geneid && exists($special_list{$geneid}) && exists($special_list{$geneid}{$id})){
			if($ar_result[1] ne "Uncertain significance" && $ar_result[2] ne "Uncertain significance")
			{
				$ar_result[0] = $ar_result[1];
			}
			if($ad_result[1] ne "Uncertain significance" && $ad_result[2] ne "Uncertain significance")
			{
				$ad_result[0] = $ad_result[1];
			}
		}
		else
		{
			if($info_ar =~ /BA1=0/){$ar_result[0] = "Benign";}
			if($info_ad =~ /BA1=0/){$ad_result[0] = "Benign";}
		}

		

		# my $note_ar = join(keys);
		# my $note_ad = ;
		print $OUT join("\t", @line[0..($#line - 1)], @ar_result, $info_ar , @ad_result, $info_ad),"\n";
		print $OUT1 join("\t", @line[0..($#line - 1)], @ar_result, @output_ar),"\n";
		print $OUT2 join("\t", @line[0..($#line - 1)], @ad_result, @output_ad), "\n";
#		return ($ar_result, $ad_result);
#		if($id eq "10_55663028_A_AT_ENST00000320301") {print Dumper(\%evidence_ad);}
#		%evidence_ad = 
	}
}

sub classfy {
    my %p_pre = @_;
    my @output = ("Uncertain significance", "Uncertain significance", "Uncertain significance");

    my @BPS = ("NA","Pathogenic","Likely pathogenic","Benign","Likely benign","Uncertain significance");
    my $PAS_out = -1; # 1:P, 2:LP
    my $BES_out = -1;# 3:B, 4:LB
    my $BPS_out = 5; # 5:Uncertain significance
    my @p = @{$p_pre{"P"}};
    my @b = @{$p_pre{"B"}};
    my $sum_p = sum @p;
    my $sum_b = sum @b;

    if($sum_p > 0){
        if($p[1] == 1){ 
            if($p[2] == 1 || $p[2] == 2){ $PAS_out = 2;}
        }
        if($p[0] >= 1){
            if($p[2] == 1){ $PAS_out = 2;}# 2:Likely pathogenic
        }
        if($p[1] == 1 && $p[3] >= 2){ $PAS_out = 2;}
        if($p[2] >= 3){$PAS_out = 2;}
        if($p[2] == 2 && $p[3] >= 2) {$PAS_out = 2;}
        if($p[2] == 1 && $p[3] >= 4) {$PAS_out = 2;}
        #-------------------------------------
        if($p[0] >= 1){
            if($p[1] >= 1) {$PAS_out = 1;} # 1:Pathogenic
            if($p[2] >= 2) {$PAS_out = 1;}
            if($p[2] == 1 && $p[3] == 1) {$PAS_out = 1;}
            if($p[3] >= 2) {$PAS_out = 1;}
        }
        if($p[0] >= 2){$PAS_out = 1;}
        if($p[1] >= 2) {$PAS_out = 1;}
        if($p[1] == 1){
            if($p[2] >= 3) {$PAS_out = 1;}
            if($p[2] == 2 && $p[3] >= 2) {$PAS_out = 1;}
            if($p[2] == 1 && $p[3] >= 4) {$PAS_out = 1;}
        }
    }
    #---------------------------------------------------------
    if($sum_b > 0){
        if($b[1] == 1 && $b[3] == 1 ) {$BES_out = 4;} #4:Likely benign
        if($b[3] >= 2 ) {$BES_out = 4;}
        #解决BP4=2的情况，将2当作3处理
        if($b[1] == 1 && $b[2] == 1 ) {$BES_out = 4;} 
        if($b[2] >= 1 && $b[3] >= 1) {$BES_out = 4;}
        #------------------------------------------------------------
        if($b[0] >= 1 || $b[1] >= 2 ) {$BES_out = 3;} #3:Benign
    }
    #-------------------------------------------------------------------------------
    if($PAS_out != -1 && $BES_out == -1) {$BPS_out = $PAS_out;}
    if($PAS_out == -1 && $BES_out != -1) {$BPS_out = $BES_out;}
    if($PAS_out == -1 && $BES_out == -1) {$BPS_out = 5;}
    if($PAS_out != -1 && $BES_out != -1) {$BPS_out = 5;}
    $output[0] = $BPS[$BPS_out];
    $output[1] = $BPS[$PAS_out];
    $output[2] = $BPS[$BES_out];
    return @output;
}
