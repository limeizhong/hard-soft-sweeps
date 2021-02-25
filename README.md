################ pi and Tajima's D for genes #############
cd input_example
python ../site2gene.pi.td.py 'example_w' 'example_gene.pos' 62

################ H12 and H2/H1 for real and simulation data #############
### real data
python ../gene.index.py example_012_pos example_gene.pos example_gene.index # get gene index first
python ../geno2h.py example_012_geno example_gene.index example.h12.stats example.h12.hap 0  # gene alone
python ../geno2h.py example_012_geno example_gene.index example.h12.stats example.h12.hap 50 # gene±50SNPs
python hap.similarity.py example_012_geno example_gene.index example.h12.sim 100 # similarity between the first two haplotypes
### simulation data
bash ../ms.sh #get simulation samples from ms simulator
python ms2h.py cul.sim.30k
python ms2h.py wild.sim.15k

################ iHH for genes #############
python ../geno2ihh.py example_012_geno example_gene.index example.gene.hap.ehh 0.05 20 

################ convert vcf to fq #############
perl ../vcf2fq.1.pl example.vcf SRR1533153 > SRR1533153.fq  # single genome
perl ../vcf2fq.2.pl example.vcf SRR1533153 SRR1533197 > SRR1533153_SRR1533197.fq # pseudo-diploid genomes
