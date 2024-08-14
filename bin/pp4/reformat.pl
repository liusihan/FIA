#!/usr/bin/perl -w

##整理文件格式：每行一个位点的所有样本-->每行一个样本的一个位点
##所有样本
die "Usage : Perl Path SourceFile DestinationFile !\n" if(@ARGV!=3);
my $dir=shift @ARGV;
my $file1="$dir/".(shift @ARGV);
my $file2="$dir/".(shift @ARGV);
open FILE1,"<$file1" or die "Can't open $file1 :$!\n";
open FILE2,">$file2" or die "Can't open $file2 :$!\n";
#print FILE2 "Sample\tCHROM\tPOS\tREF\tALT\tGT\tAD_REF\tAD_ALT\tDP\tGQ\tHom_Cntl\tHet_Cntl\tHom_Case\tHet_Case\tAF_control\tAF_case\tgnomAD_AF\tgnomAD_AF_AFR\tgnomAD_AF_AMR\tgnomAD_AF_ASJ\tgnomAD_AF_EAS\tgnomAD_AF_FIN\tgnomAD_AF_NFE\tgnomAD_AF_OTH\tgnomAD_AF_SAS\tgnomAD_Hom\tgnomAD_Hom_AFR\tgnomAD_Hom_AMR\tgnomAD_Hom_ASJ\tgnomAD_Hom_EAS\tgnomAD_Hom_FIN\tgnomAD_Hom_NFE\tgnomAD_Hom_OTH\tgnomAD_Hom_SAS\tgnomAD_Hemi\tgnomAD_Hemi_AFR\tgnomAD_Hemi_AMR\tgnomAD_Hemi_ASJ\tgnomAD_Hemi_EAS\tgnomAD_Hemi_FIN\tgnomAD_Hemi_NFE\tgnomAD_Hemi_OTH\tgnomAD_Hemi_SAS\tAC_control\tAN_control\tAC_case_all\tAN_case_all\tVEP_Consequence\tVEP_Feature\tVEP_IMPACT\tVEP_HGVSc\tVEP_HGVSp\tVEP_PICK\tVEP_Gene\tVEP_SYMBOL\n";
#更改FILE2的输出形式为：Sample	CHROM	POS	REF	ALT	VEP_Feature	VEP_SYMBOL	VEP_Gene	GT
#print FILE2 "Sample\tCHROM\tPOS\tREF\tALT\tVEP_Feature\tVEP_SYMBOL\tVEP_Gene\tGT\tFT\n";
my @title;
my @sample=();
while(my $line=<FILE1>){
	$line=~s/[\r\n]//g;
	my @line_split=split /\t/,$line;
	my $item1=join "\t",@line_split[0..3];
	my ($item2,%hash)=();
	##记录每个样本的每列信息的位置，使用数组
	my $k=0;

	if($line=~/^#/){
		for($i=0;$i<@line_split;$i++){
		#$i=7为FILE1（HLsamples_122genes.tsv）中第一个样本开始的列数
			if($line_split[$i] =~ /GT;/){
				$sample[$k]=(split /[\];:]/,$line_split[$i])[1];
				$k++;
			}
			else{
				my $colname = (split /[\];:]/,$line_split[$i])[1];
				push @title, $colname;
			}
		}
		print FILE2 join("\t","Sample", @title, "GT", "AD", "DP", "GQ","PL"), "\n";

		# print "k=$k\n";
		# print (join "\t",@sample)."\n";
	}else{
		# ##GT
		$k=0;
		##Hom_control，het_control，Hom_case，het_case
		# my ($hom_control,$het_control,$hom_case,$het_case)=(0,0,0,0);
		for($i=@title+0;$i<@line_split;$i++){
		    #$i=9为FILE1（HLsamples_122genes.tsv）中第一个样本开始的列数
		    my @info = split(/;/,$line_split[$i]);
		    my $gt = $info[0];
		    my $ad = $info[1];
		    my $dp = $info[2];
                    my $gq = $info[3];
		    my $pl = $info[4];
			if($gt ne "0/0" && $gt ne "./." && $gt ne "./0"){
				#print FILE2 "$sample[$k]\t$item1\t".(join "\t",@line_split[4..6])."\t$gt\t$ft\n";
				print FILE2 join("\t", $sample[$k], @line_split[0..$#title], $gt, $ad, $dp, $gq,$pl), "\n";
			}
			 $k++;			
		}
	}	
}
close FILE1;
close FILE2;
