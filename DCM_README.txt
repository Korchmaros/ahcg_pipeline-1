# DCM Genetic Variants

Project Goal: Find Dilated Cardiomyopathy (DCM) genetic variants in a list of six 			genes (input/1.gen_list.txt) 
Input: Bam files of 2 controls and 4 DCM patients belong to the same family     
Output: Genetic variants with reads depth and reads depth coverage plots  
     

Requirements: R version 3.1.2 and ’ggplot2’ package
$sudo apt-get install r-base-core
$R
$install.packages('ggplot2' repos='http://cran.rstudio.com', type='source')
$q()

Pipeline: Patient1
  
1. Download the BAM and BAI files.     
$ wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bam
$ wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bai      

2. Perform variant call analysis using GATK-HaplotypeCaller.
$ lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T HaplotypeCaller -R input_files/ref_genome/hg19.fa -I Patient1_RG_MD_IR_BQ.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o p1_raw_variants.vcf  

3. Recalibrate raw VCF using GATK-recalibrator.  
$ lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R input_files/ref_genome/hg19.fa -input p1_raw_variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 input_files/gatk_variantRecali/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 input_files/gatk_variantRecali/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 input_files/gatk_variantRecali/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 input_files/dbsnp/dbsnp_138.hg19.vcf -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile p1_recalibrate_SNP.recal -tranchesFile p1_recalibrate_SNP.tranches

$ lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R input_files/ref_genome/hg19.fa -input p1_raw_variants.vcf -mode SNP --ts_filter_level 99.0 -recalFile p1_recalibrate_SNP.recal -tranchesFile p1_recalibrate_SNP.tranches -o p1_recalibrated_snps_raw_indels.vcf

4. Create the BED file of the coordinates for each of the DCM genes. The output is temp/exom_list.bed.     
$ perl scripts/grep_gen_list.pl

5. Extract the variants associated to each of the DCM genes.    
$ bedtools intersect -wa -header -a p1_recalibrated_snps_raw_indels.vcf -b scripts/temp/exom_list.bed > patient1_dcm_final.vcf

6. Calculate the reads depth information for each of the DCM genes and create the reads depth plot.
   6.1 Extract the alignments using samtools.    
$ samtools view -L scripts/temp/exom_list.bed Patient1_RG_MD_IR_BQ.bam -b > patient1.dcm.bam
   6.2 Compute and summarize the depth for each of the DCM genes.    
$ bedtools genomecov -ibam patient1.dcm.bam -bga > patient1.dcm.bed
   6.3 Intersect patient1.dcm.bed and exom_list.bed. 
$ bedtools intersect -split -a patient1.dcm.bed -b scripts/temp/exom_list.bed -bed > patient1.final.bed
   6.4 Calculate the coverage depth.
$ python scripts/depth_for_each_site.py
$ perl scripts/label.pl scripts/temp/patient1.final.full.txt scripts/txt/p1.final.txt

7. For each of DCM genes, find the depth coverage cutoffs and plot the depth coverage distribution of each exon.
$ R CMD BATCH '--args scripts/txt/p1.final.txt' scripts/plots.R
$ cp boxplot.exons.LMNA.jpg SNP_report/boxplots/P1.boxplot.exons.LMNA.jpg
$ cp boxplot.exons.MYBPC3.jpg SNP_report/boxplots/P1.boxplot.exons.MYBPC3.jpg
$ cp boxplot.exons.MYH6.jpg SNP_report/boxplots/P1.boxplot.exons.MYH6.jpg
$ cp boxplot.exons.MYH7.jpg SNP_report/boxplots/P1.boxplot.exons.MYH7.jpg
$ cp boxplot.exons.SCNSA.jpg SNP_report/boxplots/P1.boxplot.exons.SCNSA.jpg
$ cp boxplot.exons.TNNT2.jpg SNP_report/boxplots/P1.boxplot.exons.TNNT2.jpg
$ rm boxplot.exons.LMNA.jpg boxplot.exons.MYBPC3.jpg boxplot.exons.MYH6.jpg boxplot.exons.MYH7.jpg boxplot.exons.SCNSA.jpg boxplot.exons.TNNT2.jpg
$ cp cutoff.MYH7.csv SNP_report/cutoffs/P1.cutoff.MYH7.csv
$ cp cutoff.MYH6.csv SNP_report/cutoffs/P1.cutoff.MYH6.csv
$ cp cutoff.LMNA.csv SNP_report/cutoffs/P1.cutoff.LMNA.csv
$ cp cutoff.TNNT2.csv SNP_report/cutoffs/P1.cutoff.TNNT2.csv
$ cp cutoff.SCNSA.csv SNP_report/cutoffs/P1.cutoff.SCNSA.csv
$ cp cutoff.MYBPC3.csv SNP_report/cutoffs/P1.cutoff.MYBPC3.csv
$ rm cutoff.MYH7.csv cutoff.MYH6.csv cutoff.LMNA.csv cutoff.TNNT2.csv cutoff.SCNSA.csv cutoff.MYBPC3.csv

8. Find the variants that have a depth above the cutoff.
$ cp cutoff.final.csv scripts/cutoff/P1cutoff.final.csv
$ rm cutoff.final.csv
$ bedtools intersect -wa -a p1_recalibrated_snps_raw_indels.vcf -b scripts/temp/exom_list.bed > scripts/vcf/P1_afterinter.vcf
$ python scripts/min.py
$ cp scripts/vcf/*new.vcf SNP_report/

9. Update ‘SNP_report’ folder.
$ rm scripts/vcf/*new.vcf
$ rm plots.Rout Patient1_RG_MD_IR_BQ.bam p1_raw_variants.vcf patient1.dcm.bam p1_recalibrate_SNP.recal Patient1_RG_MD_IR_BQ.bai p1_recalibrated_snps_raw_indels.vcf.idx p1_raw_variants.vcf.idx patient1.dcm.bed patient1.dcm.bed patient1_dcm_final.vcf p1_recalibrate_SNP.tranches.pdf p1_recalibrate_SNP.tranches.pdf p1_recalibrate_SNP.tranches p1_recalibrate_SNP.recal.idx p1_recalibrated_snps_raw_indels.vcf

Repeat those 9 steps for all the patients and the controls.

Find the variants that are more likely to be deleterious.
$ perl scripts/fliter.pl

Output Folder(SNP_Report):
1)The final output is SNP_report.txt
2)Boxplot folder contains 36 boxplots(one for each gene and for each person)
3)Cutoff folder contains 36 tables with median, 90th quantile, and cutoffs of the depth (one for each gene and for each person)   

Methodology:
For each gene and for each exon the cutoff is the lower whisker of the depth distribution.
lower whisker = Q_1 - 1.5 (Q3-Q1), Q3=75-th quantile Q1=25-th quantile

In p1_recalibrated_snps_raw_indels.vcf the variants are classified according to
 
if cutoff-depth < 0, “Mutation untrusted” 
if 0< cutoff-depth < 5, “Ambiguous sequencing result” 
if cutoff-depth > 4, “Mutation confirmed” 

In scripts/fliter.pl the variants more likely to cause DCM are selected according to these scores. 

Patients mutations
“Mutation untrusted” score=1
“Ambiguous sequencing result” score=0.6
“Mutation confirmed” score=0

Controls mutations
“Mutation confirmed” presented in none of the controls, score=-1
“Mutation confirmed” presented in one of the controls, score=1
“Mutation confirmed” presented in both controls, not considered

Authors:
ANNACHIARA KORCHMAROS 
CHEN GUO
MENGNAN ZHANG
PEIJUE ZHANG


