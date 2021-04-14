############### PART1 DEMOGRAPHIC ANALYSES ################
# 1 download SRA raw reads
prefetch srr_accession
fastq-dump -gzip --split-3 srr_accession

# 2 creat index for the reference genome 
samtools faidx reference_genome.fa
bwa index reference_genome.fa

# 3 align to reference genome and sort
bwa mem -t 48 -M -R "@RG\tID:srr_accession\tSM:srr_accession\tLB:WGS\tPL:Illumina" reference_genome.fa srr_accession_1.fastq.gz srr_accession_2.fastq.gz | samtools sort -@ 16 -O BAM -o srr_accession.bam

# 4 mark the duplicate regions and creat index 
gatk MarkDuplicates --java-options '-Xmx30G' -I srr_accession.bam -O srr_accession.markdup.bam -M srr_accession.metrics 
samtools index srr_accession.markdup.bam

# 5 generate gvcf files
gatk HaplotypeCaller --java-options '-Xmx20G -XX:ParallelGCThreads=8' --emit-ref-confidence GVCF -R reference_genome.fa -I srr_accession.markdup.bam -L chromosome -O prefix0.g.vcf.gz
gatk CombineGVCFs -R reference_genome.fa gvcfs_for_all_samples -O prefix1.g.vcf.gz

# 6 SNPs calling and filtering
gatk GenotypeGVCFs -R reference_genome.fa -V prefix1.g.vcf.gz --all-sites -O prefix1.all.vcf.gz
gatk GenotypeGVCFs -R reference_genome.fa -V prefix1.g.vcf.gz -O prefix1.var.vcf.gz
bcftools filter -O z -o prefix1.all.filter.vcf.gz -e 'QD< 2.0 || FS> 60.0 || MQ< 40.0' prefix1.all.vcf.gz
bcftools filter -O z -o prefix1.var.filter.vcf.gz -e 'QD< 2.0 || FS> 60.0 || MQ< 40.0' prefix1.var.vcf.gz

# 7 Convert vcf files to fastq files that PSMC required
perl vcf2fq.2.pl prefix1.all.filter.vcf.gz srr_accession_1 srr_accession_2 > srr1.srr2.fq　

# 8 PSMC analyses
fq2psmcfa srr1.srr2.fq > srr1.srr2.psmcfa
psmc -p "4+25*2+4+6" -o srr1.srr2.psmc srr1.srr2.psmcfa
psmc_plot.pl -x 500 -u 6e-09 -g 1 -p srr1.srr2.plot srr1.srr2.psmc 

# 9 SMC++ analyses,split model
smcpp='docker run --rm -v your/path/way/:/mnt terhorst/smcpp:latest'
$smcpp vcf2smc prefix1.var.filter.vcf.gz data_w/c.prefix1.smc.gz chromosome Wild:sample1,sample2,...
$smcpp vcf2smc prefix1.var.filter.vcf.gz data_c/c.prefix1.smc.gz chromosome cultivated:sample1,sample2,...
$smcpp estimate -o pop_w/ 6e-9 data_w/*.smc.gz
$smcpp estimate -o pop_c/ 6e-9 data_c/*.smc.gz
cp data_w/*.smc.gz data_wc/
cp data_c/*.smc.gz data_wc/
$smcpp split -o split/ pop_w/model.final.json pop_c/model.final.json data_wc/*.smc.gz
$smcpp plot -g 1 -c wild_cultivated.pdf split/model.final.json


############### PART2 SELECTION ANALYSES ################ 
#calculate pi and Tajima's D for genes
cd input_example
python ../site2gene.pi.td.py 'example_w' 'example_gene.pos' 62

#calculate H12 and H2/H1 for real genome data
python ../gene.index.py example_012_pos example_gene.pos example_gene.index 
python ../geno2h.py example_012_geno example_gene.index example.h12.stats example.h12.hap 0  # gene alone
python ../geno2h.py example_012_geno example_gene.index example.h12.stats example.h12.hap 50 # gene±50SNPs
python hap.similarity.py example_012_geno example_gene.index example.h12.sim 100 # similarity between the first two haplotypes

#calculate H12 and H2/H1 for simulation data
bash ../ms.sh #get simulation samples from ms simulator
python ms2h.py cul.sim.30k
python ms2h.py wild.sim.15k

#calculate iHH for genes
python ../geno2ihh.py example_012_geno example_gene.index example.gene.hap.ehh 0.05 20

#calculate Fst between two populations for each SNP
vcftools --vcf input.vcf --weir-fst-pop pop1 --weir-fst-pop pop2 --out pop1_pop2
