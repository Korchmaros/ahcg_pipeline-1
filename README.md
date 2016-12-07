# DCM Genetic Variants

Project Goal: Find Dilated Cardiomyopathy (DCM) related genetic variants      
Input: Bam files of 2 controls and 4 DCM patients      
Output: Genetic variants with reads depth and reads depth coverage plot  

Procedures:     
1) Perform variant call analysis and variant recalibration     
2) Read depth calculation on DCM-related gene panel       
3) Read depth coverage plot for DCM-related gene panel       

Requirements: R version 3.1.2 and ’ggplot2’ package
$sudo apt-get install r-base-core
$R
$install.packages('ggplot2' repos='http://cran.rstudio.com', type='source')
$q()

Example: Patient1
  
1. Download the BAM and BAI files.     
$ wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bam
$ wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bai      

2. Perform variant call analysis using GATK-HaplotypeCaller.
$ lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T HaplotypeCaller -R input_files/ref_genome/hg19.fa -I Patient1_RG_MD_IR_BQ.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o p1_raw_variants.vcf  

3. Recalibrate raw VCF using GATK-recalibrator.  
$ lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R input_files/ref_genome/hg19.fa -input p1_raw_variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 input_files/gatk_variantRecali/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 input_files/gatk_variantRecali/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 input_files/gatk_variantRecali/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 input_files/dbsnp/dbsnp_138.hg19.vcf -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile p1_recalibrate_SNP.recal -tranchesFile p1_recalibrate_SNP.tranches

$ lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R input_files/ref_genome/hg19.fa -input p1_raw_variants.vcf -mode SNP --ts_filter_level 99.0 -recalFile p1_recalibrate_SNP.recal -tranchesFile p1_recalibrate_SNP.tranches -o p1_recalibrated_snps_raw_indels.vcf

4. Create the BED file of the coordinates for each of the DCM genes. It creates in the folder temp the file exom_list.bed    
$ perl scripts/grep_gen_list.pl

5. Extract the variants associated to each of the DCM genes.    
$ bedtools intersect -wa -header -a p1_recalibrated_snps_raw_indels.vcf -b scripts/temp/exom_list.bed > patient1_dcm_final.vcf    

6. Calculate the reads depth information for of each of the DCM genes and create the reads depth plot.
   6.1 Extract the alignments using samtools.    
$ samtools view -L scripts/temp/exom_list.bed Patient1_RG_MD_IR_BQ.bam -b > patient1.dcm.bam
   6.2 Compute and summarize the depth for each of the DCM genes.    
$ bedtools genomecov -ibam patient1.dcm.bam -bga > patient1.dcm.bed
   6.3 Intersect patient1.dcm.bed and exom_list.bed and calculate the reads’ depth. 
$ bedtools intersect -split -a patient1.dcm.bed -b scripts/temp/exom_list.bed -bed > patient1.final.bed
$ python scripts/depth_for_each_site.py
$ cp patient1.final.full.labeled.txt plotting/p1.final.txt
$ rm patient1.final.full.labeled.txt

7. For each of DCM genes, find the depth coverage cutoffs and plot the depth coverage distribution of each exon.
$ R CMD BATCH '--args scripts/txt/p1.final.txt' scripts/plots.R
$ cp boxplot.exons.LMNA.jpg SNP_report/boxplots
$ cp boxplot.exons.MYBPC3.jpg SNP_report/boxplots
$ cp boxplot.exons.MYH6.jpg SNP_report/boxplots
$ cp boxplot.exons.MYH7.jpg SNP_report/boxplots
$ cp boxplot.exons.SCNSA.jpg SNP_report/boxplots
$ cp boxplot.exons.TNNT2.jpg SNP_report/boxplots
$ rm boxplot.exons.LMNA.jpg boxplot.exons.MYBPC3.jpg boxplot.exons.MYH6.jpg boxplot.exons.MYH7.jpg boxplot.exons.SCNSA.jpg boxplot.exons.TNNT2.jpg
$ cp cutoff.MYH7.csv SNP_report/cutoffs
$ cp cutoff.MYH6.csv SNP_report/cutoffs
$ cp cutoff.LMNA.csv SNP_report/cutoffs
$ cp cutoff.TNNT2.csv SNP_report/cutoffs
$ cp cutoff.SCNSA.csv SNP_report/cutoffs
$ cp cutoff.MYBPC3.csv SNP_report/cutoffs
$ rm cutoff.MYH7.csv cutoff.MYH6.csv cutoff.LMNA.csv cutoff.TNNT2.csv cutoff.SCNSA.csv cutoff.MYBPC3.csv

8. Final analysis
8.1 Eliminate a variant if its depth is below the cutoff, leave a note if the difference between the depth and the cutoff is below 5.
$ cp P1cutoff.final.csv scripts/P1cutoff.final.csv
$ rm P1cutoff.final.csv
$ bedtools intersect -wa -header -a p1_recalibrated_snps_raw_indels.vcf -b scripts/temp/exom_list.bed > scripts/vcf/P1_afterinter.vcf
$ python scripts/min.py
$ cp scripts/vcf/*new.vcf SNP_report/
8.2 Update ‘SNP_report’ output folder and remove unnecessary files
$ rm scripts/vcf/*new.vcf
$ rm plots.Rout Patient1_RG_MD_IR_BQ.bam p1_raw_variants.vcf patient1.dcm.bam p1_recalibrate_SNP.recal Patient1_RG_MD_IR_BQ.bai p1_recalibrated_snps_raw_indels.vcf.idx p1_raw_variants.vcf.idx patient1.dcm.bed patient1.dcm.bed patient1_dcm_final.vcf p1_recalibrate_SNP.tranches.pdf p1_recalibrate_SNP.tranches.pdf p1_recalibrate_SNP.tranches p1_recalibrate_SNP.recal.idx p1_recalibrated_snps_raw_indels.vcf


The outputs are in ‘SNP_report’ folder


