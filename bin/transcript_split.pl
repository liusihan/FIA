#!/usr/bin/perl -w

die "Usage: Perl dir :$!" if(@ARGV!=1);
my $dir=shift @ARGV;
my $file1="$dir/SNV_allsites_v1.rmsk.tsv";
my $file2="$dir/SNV_allsites_v2.rmsk.tsv";
open FILE1,"<$file1" or die "Can't open $file1 :$!\n";
open FILE2,">$file2" or die "Can't open $file2 :$!\n";
print FILE2 "CHROM\tPOS\tREF\tALT\tGNOMAD_AF_SAS_WES\tgnomAD_grpmax_WES\tgnomAD_grpmax_WGS\tgnomAD_AF_EAS_WES\tgnomAD_AF_ASJ_WES\tgnomAD_AF_MID_WES\tgnomAD_AF_FIN_WES\tgnomAD_AF_NFE_WES\tgnomAD_AF_NFE_WES\tDBNSFP_INTERPRO_DOMAIN\tVEP_MaxEntScan_diff\tVEP_Consequence\tVEP_DOMAINS\tVEP_Feature\tVEP_IMPACT\tVEP_HGVSc\tVEP_HGVSp\tVEP_PICK\tVEP_Gene\tgnomAD_AF_AMI_WES\tclinvar_CLNSIG\tDVD_PATHOGENICITY\tDBNSFP_REVEL_SCORE\trmsk_repName\n";

my @VepIndex=();
my $col=0;
while(my $line=<FILE1>){
	$line=~s/[\r\n]//g;
	my @line_split=split /\t/,$line;
	if($line=~/^#/){
		my $k=0;
		foreach my $i(0..$#line_split){
			$line_split[$i]=~s/[\[\]0-9]+//g;
			if($line_split[$i]=~/VEP_/){
				$VepIndex[$k]=$i ;
				$k++;
			}
		}
		$col=$#line_split;
	}else{
		foreach my $i(0..$col){
			$line_split[$i]="" if(not defined $line_split[$i]);
			# my $item=(join "\t",@line_split[0..14,24..27]);
		}
		my $num=()=$line_split[$VepIndex[0]]=~/\|/g;
		foreach my $i(0..$num){
			print FILE2 (join "\t",@line_split[0..13])."\t";
			foreach my $k(0..$#VepIndex){
				my @array=split /\|/,$line_split[$VepIndex[$k]];
				$array[$i]="" if(not defined $array[$i]);
				print FILE2 $array[$i]."\t";
			}
			print FILE2 (join "\t",@line_split[23..27])."\n";
		}
	}
}
close FILE1;
close FILE2;
