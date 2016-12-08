How to calculate read depth based on alignment file.    
(All the files in the folder "cancer_gene_panel")    

1. Extract cancer gene chromosome coordinates for the given gene list  
1.1 Extract gene chromosome coordinates by matching transcript ID.  
`$ python get_coordinates_from_gene.py BRC_OC_gene_list.txt hg19_refGene.txt > brc_oc_gene_list_bed.txt`    

1.2. Add 20bp to both ends of chromosome coordinates for each gene.    
`$ python 20bp_flank_both_ends.py brc_oc_gene_list_bed.txt > brc_oc_gene_list_bed_add_20_flank.txt`    
1.3 Add exon information
`$ python add_exon_to_nm.py brc_oc_gene_list_bed_add_20_flank.txt > brc_oc_gene_list_bed_add_20_flank_add_exon.txt`  

2. Use bedtools to cut the regions from vcf file.  
$ bedtools intersect -wa -a <vcf file> -b <bed file>  
`$ bedtools intersect -wa -a recalibrated_snps_raw_indels.vcf -b brc_oc_gene_list_bed_add_20_flank.txt > brc_oc_gene_list.vcf`


3. Calculate the coverage.
3.1. Extract brca1 alignments

http://vannberg.biology.gatech.edu/data/ahcg2016/NA12878_hg19_final.bam 
http://vannberg.biology.gatech.edu/data/ahcg2016/NA12878_hg19_final.bai
http://vannberg.biology.gatech.edu/data/ahcg2016/vcf/NA12878_variants.vcf

`$ samtools view -L brca1.bed NA12878_hg19_final.bam -b > NA12878.brc_oc_gene_list.bam`
      Note: -L: only output alignments overlapping in the input bed file
      Note: -b: output alignments in the bam format

3.2. Computes and summarize coverage for brca1
`$ bedtools genomecov -ibam NA12878.brc_oc_gene_list.bam -bga > NA12878.brc_oc_gene_list.bga.bed`
      Note: -ibam BAM file as input for coverage.
      Note: -bga Reporting genome coverage for all positions in BEDGRAPH format.
      Note: This is BedGraph format:
                chrom chromStart chromEnd dataValue

3.3. Intersection between two bed files.
`$ bedtools intersect -split -a NA12878.brc_oc_gene_list.bga.bed -b brc_oc_gene_list_bed_add_20_flank.txt -bed > brc_oc_gene_list.final.bed`
      Note: when using the -split option, only the exon overlaps are reported  
      
3.4. Add exon information
`$ python add_exon_to_depth.py brc_oc_gene_list_bed_add_20_flank_add_exon.txt brc_oc_gene_list.final.bed > brc_oc_gene_list.exon.final.bed`

3.5. Use the ReadDepthCalculation.py to get the read depth result.
`$ python read_depth_calculation_w_exon.py brc_oc_gene_list.exon.final.bed > brc_oc_depth_exon.txt`


4. Pathogenic annotation for the variants.
4.1. Download the BRCA variant pathogenic annotation file.
$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA1_brca_exchange_variants.csv
$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA2_brca_exchange_variants.csv

4.2. Use the PathgenicAnnotation.py script to match the variants between vcf file and pathogenic annotation file.
$ python pathgenic_annotation_add_dp.py recalibrated_snps_raw_indels.vcf BRCA2_brca_exchange_variants.csv > BRCA2_variants_annotation_add_dp.vcf


5. Integrate the reads coverage information with the pathogenic information for variants.(Need to specify the cutoff as sys.argv[3])
$ python report_depth_pathogenic.py BRCA2_variants_annotation_add_dp.vcf brc_oc_depth_exon.txt 10 > brc_oc_depth_exon_confidence.txt


6. Draw the coverage plot.
`$ Rscript draw_depth.R brc_oc_depth_exon.txt`  
