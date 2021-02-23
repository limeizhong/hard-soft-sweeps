#!/usr/bin/perl -w
use strict;
use constant MIN_DP => 10;
use constant MAX_DP => 100;

my %het = (AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y',
             GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K',
			 AA=>'A', TT=>'T', GG=>'G',CC=>'C');
			 
my $infile = shift;
my $id1=shift;
my $id2=shift;
open(IN, "gzip -c -d $infile|") or die "Cannot open $infile: $!";
my @index=();
my $last_pos=1;
my $last_chr=undef;
my ($seq, $qual)=();
while(<IN>){
	chomp;
	if(/^#{2,}/){
		next;
	}elsif(/^#CHR/){
		my @headers=split/\s+/;
		#print join("\t", @headers[0...8]),$id1.'_'.$id2,"\n";
		splice(@headers, 0, 9);
		foreach (0..$#headers){
			if($headers[$_] eq $id1 or $headers[$_] eq $id2){
				push @index, $_;
			}
		}  
		# print STDERR join(', ',@index), " were chosen\n";
	}else{
		my @inds=split/\s+/;
		my ($chr, $pos, $ref, $alt, $filter, $info,$format)=@inds[0,1,3,4,6,7,8];
		$last_chr = $chr unless(defined $last_chr);
		if($chr ne $last_chr){ # a new chr
			$last_chr = $chr;
			$last_pos = 1;
			print "\@$chr\n";
			print_seq($seq);
			print "+\n";
			print_seq($qual);
		}
		
		if($pos - $last_pos >1){
			$seq.='n' x ($pos - $last_pos - 1);
			$qual.= '!' x ($pos - $last_pos - 1);
		}
		
		$last_pos = $pos;
		
		my ($MQ) = $info=~/MQ=(\d+)/;
		$MQ ||=0; 
		$MQ = 126 if($alt eq '.');
		splice(@inds,0,9);
		
		if($filter eq 'FILTERED'){
			$seq.= 'n';
			$qual.='!';
		}else{
			my $gt1_str = $inds[$index[0]];
			my $gt2_str = $inds[$index[1]];
			my $allele1=(split/\/|\|/, (split/:/,$gt1_str)[0])[0];
			my $allele2=(split/\/|\|/, (split/:/,$gt2_str)[0])[0]; #两个不同的个体
			my $ad_str1 = (split/:/,$gt1_str)[1]; #REF和支持ALT的测序深度。
			my $ad_str2 = (split/:/,$gt2_str)[1];
			
			if($allele1 ne '.' && $allele2 ne '.' && defined $ad_str1 && defined $ad_str2 && $ad_str1=~/\d/ && $ad_str2=~/\d/){
				my ($base1, $base2, $dp1, $dp2);
				
				my @alleles = ($ref, (split/,/,$alt));
				$base1 = $alleles[$allele1];
				$base2 = $alleles[$allele2];
				# print $base1,$base2,", ";
				if(length($base1)>1 || length($base2)>1 || $base1 eq '*' || $base2 eq '*'){ #indel
					$seq.='n';
					$qual.='!';
				}else{
					$dp1 = (split/,/,$ad_str1)[$allele1];
					$dp2 = (split/,/,$ad_str2)[$allele2];
					$base1 = uc($base1);
					$base2 = uc($base2);
					my $genotype = $het{$base1.$base2};
					$genotype ||='N';
					if($dp1 + $dp2 < MIN_DP || $dp1 + $dp2 > MAX_DP || $filter eq 'INDEL_FLANK' ){
						$genotype=lc($genotype);
					}
					$seq.=$genotype;
					my $q =int($MQ + 33 + .499);
					$q = chr($q <= 126? $q : 126);
					
					$qual.=$q;
				}
			}else{ #  no genotype
				$seq.='n';
				$qual.='!';
			}
		}
	}
}
close(IN);
# last chr
print "\@$last_chr\n";
print_seq($seq);
print "+\n";
print_seq($qual);	
	
sub print_seq{
	my ($s) = @_;
	my $l = length($s);
  for (my $i = 0; $i < $l; $i += 60) {
    print substr($s, $i, 60), "\n";
  }
}
	