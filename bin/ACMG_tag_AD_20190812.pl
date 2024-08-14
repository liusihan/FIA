#!/usr/bin/perl -w

die "Usage: Perl dir PathoFile DVPred:$!" if(@ARGV!=2);
my $dir=shift @ARGV;
chdir $dir or die "Can't chdir to $dir :$!\n";

my $file1=shift @ARGV;

open FILE1,"<$file1" or die "Can't open $file1 :$!\n";
my %hash1=();
while(my $line=<FILE1>){
	$line=~s/[\r\n]//g;
	my @line_split=split /\t/,$line;
	my $item=join ":",@line_split[0,1,2,3];
	$hash1{$item}=$line_split[4];
}
close FILE1;


$file1="SNV_allsites_v2.rmsk.tsv";
$file2="SNV_allsites_ACMG_AD.tsv";
open FILE1,"<$file1" or die "Can't open $file1 :$!\n";
open FILE2,">$file2" or die "Can't open $file2 :$!\n";
print FILE2 "chrom\tpos\tref\talt\tgnomAD_AF_SAS_WES\tgnomAD_AF_SAS_WES\tgnomAD_grpmax_WES\tgnomAD_AF_grpmax_WGS\tgnomAD_AF_EAS_WES\tgnomAD_AF_ASJ_WES\tgnomAD_AF_MID_WES\tgnomAD_AF_FIN_WES\tgnomAD_AF_NFE_WES\tgnomAD_AF_AFR_WES\tDBNSFP_INTERPRO_DOMAIN\tVEP_MaxEntScan_diff\tVEP_CONSEQUENCE\tVEP_DOMAINS\ttranscript\tGNOMAD_AF_AMI_WES\tVEP_HGVSC\tVEP_HGVSP\tVEP_PICK\tVEP_GENE\tVEP_SYMBOL\tCLINVAR_CLNSIG\tDVD_PATHOGENICITY\tDBNSFP_REVEL_SCORE\trmsk_repName\tPP5\tPM2\tBA1\tBS1\tPM4\treplace1\treplace2\tBP3\tPP7\tBS2\n";
my ($pp5,$pm2_2,$pm2_3,$ba1,$bs1, $bs2,$pm4,$pp3,$bp4,$pp7,$bp3)=();
while(my $line=<FILE1>){
	$line=~s/[\r\n]//g;
	next if($line=~/CHROM/);
	my @line_split=split /\t/,$line;
	next if($line_split[3] eq "*");
	foreach my $i(0..37){
		$line_split[$i]="" if(not defined $line_split[$i]);
	}
	my $item=join ":",@line_split[0,1,2,3];
	
	##PP5=0
	if(exists $hash1{$item}){
#		if($hash1{$item}=~/P/){
			$line_split[28]=0 ;
			$pp5++;
	}
#		}
#		$line_split[3]="$line_split[3]\t$ha";
#	}else{
		$line_split[3]="$line_split[3]\t";
#	}
	my $maxaf = 0;
	my @afcols=(4..11);
	push @afcols,23;
	foreach my $i(@afcols){
		if($line_split[$i] =~ /\d+/ && $line_split[$i] > $maxaf){
			$maxaf = $line_split[$i];
		}
	}

	if($maxaf >=0.001){
			$line_split[30]=0;
			$ba1++;
		}
	elsif($maxaf >=0.0002 && $maxaf <0.001){
		$line_split[31]=1;
		$bs1++;
		}
	elsif($maxaf <=0.00002){
		$line_split[29]=3;
		$pm2_3++;
		}
	
	##PM4=2
	if(!($line_split[27]=~/\w+/) && $line_split[15]=~/inframe/){
		$line_split[32]=2;
		$pm4++;
	}

	if($line_split[27]=~/\w+/ && $line_split[15]=~/inframe/ && !($line_split[16] =~/\w+/) ){
		$line_split[35]=3;
		$bp3++;
	}
	
	if($line_split[26] ne "."){
		if($line_split[26] >= 0.7){
			$line_split[33]=3;
			$pp3++;
		}
		#BP4=3
		if($line_split[26] <= 0.15){
			$line_split[34]=3;
			$bp4++;
		}
	}
#	}else{
		#$line_split[28]="\t$line_split[28]";
#	}
	
	##BS2=1
	if($line_split[12]=~/\d+/ && $line_split[12] >=1){
		$line_split[37]=1;
		$bs2++;
	}
	## PP7=3
	if($line_split[14]=~/\d+/ && $line_split[14]>=6){
		$line_split[33]=3;
		$pp3++;
	}

	print FILE2 (join "\t",@line_split)."\n";
}
close FILE1;
close FILE2;
# print "---------AD--------\n";
# printf ("PP5=0: %d\n",$pp5);
# printf ("PM2=2: %d\n",$pm2_2);
# printf ("BA1=0: %d\n",$ba1);
# printf ("BS1=1: %d\n",$bs1);
# printf ("PM4=2: %d\n",$pm4);
# printf ("BP3=3: %d\n",$bp3);
# printf ("PP3=2: %d\n",$pp3);
# printf ("BP4=3: %d\n",$bp4);
# printf ("BS2=1: %d\n",$bs2);
# printf ("DVPred>0.5: %d\n",$DVPred);
# print "-------------------\n";
