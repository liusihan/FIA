use warnings;
use strict;
use Data::Dumper;
use List::Util qw/sum/;
use File::Basename;
use Storable qw(dclone);

die "perl $0 <sample_file> <outdir> \n" unless(@ARGV == 2);
my $samplefile = $ARGV[0];
my $outdir = $ARGV[1];
# my $samplelist = $ARGV[2];
# # my $outdir = $ARGV[3];
# open(FILE, $samplelist) or die $!;
# my $content = do {local $/=undef;<FILE>;};
# my @samplelist = split(/\n/, $content);


#foreach my $sample(@samplelist){
	my $sample_prefix = basename($samplefile, ".classfy.tsv");
	my %p_allele_ad;
	my %p_allele_ar;
	#print $sample_prefix,"\n";
	open(my $H, $samplefile);
	countP($H);
	close $H;
	#print Dumper(\%p_allele_ar),"\n";
	open(my $H1, $samplefile);
	open(my $OUT, ">$outdir/${sample_prefix}.classfy.Pallele.tsv");
	Pallele($H1, $OUT);
	close $H1;
	close $OUT;
#}


sub countP {
	my ($H) = @_;
	my @indx = (0,0,0);

	while(<$H>){
		chomp;
		my @line = split(/\t/, $_);
		if($line[0] eq "chrom"){
			$indx[0] = (grep {$line[$_] eq "AR_classify"} 0..$#line)[0];
			$indx[1] = (grep {$line[$_] eq "AD_classify"} 0..$#line)[0];
			$indx[2] = (grep {$line[$_] eq "gt"} 0..$#line)[0];
			next;
		}
		my $id = join("_", @line[0..4]);
#		print $id,"\n";

		if($line[$indx[0]] eq "Pathogenic" || $line[$indx[0]] eq "Likely pathogenic"){
			if($line[$indx[2]] eq "1/1" || $line[$indx[2]] eq "./1"){
				$p_allele_ar{$line[4]}+=2;
			}
			else{
				$p_allele_ar{$line[4]}++;
			}
		}
		if($line[$indx[1]] eq "Pathogenic" || $line[$indx[1]] eq "Likely pathogenic"){
			if($line[$indx[2]] eq "1/1" || $line[$indx[2]] eq "./1"){
				$p_allele_ad{$line[4]}+=2;
			}
			else{
				$p_allele_ad{$line[4]}++;
			}
		}
	}
}



sub Pallele {
	my ($H1, $OUT1) = @_;
	while(<$H1>){
		chomp;
		my @line = split(/\t/, $_);
		if($line[0] eq "chrom"){print $OUT1 join("\t", "Sample", $_, "P_number_AR", "P_number_AD"),"\n"; next;}
		my $id = join("_", @line[0..4]);
#		print $id,"\n";
		my @p_allele = (0,0);
		if(exists($p_allele_ar{$line[4]})){
			$p_allele[0] = $p_allele_ar{$line[4]};
		}
		if(exists($p_allele_ad{$line[4]})){
			$p_allele[1] = $p_allele_ad{$line[4]};
		}
		print $OUT1 join("\t", $sample_prefix, @line, @p_allele),"\n";
	}
}
