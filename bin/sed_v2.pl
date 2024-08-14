#!/usr/bin/perl -w

##只输出以下列：chrom,pos,ref,alt,transcript,PP5,PM2,BA1,BS1,PM4,PP3,BP4,BP3,PP7,BS2
# my $dir="/public/home/swgenetics/task/20180607/ACMG/output";
die "Usage: Perl dir :$!" if(@ARGV!=2);
# my $dir_ar=shift @ARGV;
# my $dir_ad=shift @ARGV;
my @dir = @ARGV;
#chdir $dir or die "Can't chdir to $dir :$!\n";
my @array=("SNV_allsites_ACMG_AD.tsv", "SNV_allsites_ACMG_AR.tsv");

foreach my $i(0..$#array){
	my $file1 = $array[$i];
	open FILE1,"<$file1" or die "Can't open $file1 :$!\n";
	my $file2=$array[$i];
	$file2=~s/tsv/stdout.tsv/;
	open FILE2,">$dir[$i]/$file2" or die "Can't open $file2 :$!\n";
	while(my $line=<FILE1>){
		$line=~s/[\r\n]//g;
		my @line_split=split /\t/,$line;
		my $num=()=$line=~/\t/g;
		# $line_split[18]="" if(not defined $line_split[18]);
		my $item="$line_split[0]\t$line_split[1]\t$line_split[2]\t$line_split[3]\t$line_split[18]";
		foreach my $i(29..38){
			$line_split[$i]="NA" if(not defined $line_split[$i] or !($line_split[$i] =~/\d/));			
		}		
		print FILE2 $item."\t".(join "\t",@line_split[29..38])."\n";
	}
	close FILE1;
	close FILE2;
}
