GTFtools
====




## Description

GTFtools provides a set of functions to compute or extract various features of gene models as described in the table below. Note that GTFtools can be applied to not only human but also non-human gene models like the lab mouse.

| **Options** | **Functions and example use**                          | **Notes** |
|-------------|--------------------------------------------------------|-----------|
| -h | Help information.<BR/>Example use:<BR/>gtftools -h |           |
| -m | For each gene, calculate merged exons by merging exons of all splice isoforms from the same gene. The output is merged exons in bed format.<BR/>Example use:<BR/>gtftools -m merged_exons.bed demo.gtf |Used to calculate nonverlapping exonic length of genes with multiple splice isoforms.|
| -d | Calculate independent introns which is defined as introns (or part of introns) that do not overlap with any exons of any genes in the genome2. It is calculated by subtracting merged exons from genes. The output is in bed format.<BR/>Example use:<BR/>gtftools -d independent_introns.bed demo.gtf |Used in intron retention detection.|
| -l | Calculate gene lengths. Since a gene may have multiple isoforms, there are multiple ways to calculate gene length based on literature. Three simple ways are considering the mean, median and maximum of the lengths of isoforms as the length of the gene. A fourth way is to calculate the length of merged exons of all isoforms (i.e. non-overlapping exonic length). So, in total, four different types of gene lengths(the mean, median and max of lengths of isoforms of agene, and the length of merged exons of isoforms of a gene) are provided. format.<BR/>Example use:<BR/>gtftools -l gene_length.txt demo.gtf |Needed for<BR/> e.g.<BR/> calculating FPKM in RNA-seq data analysis, where gene length is required.|
| -r | Calculate transcript isoform lengths.<BR/>Example use:<BR/>gtftools -r isoform_length.txt demo.gtf||
| -g | Output gene coordination and ID mappings in bed format.<BR/>Example use:<BR/>gtftools -g genes.bed demo.gtf||
| -p | An input file containing a list of SNPs with at least three columns, with the first being chromosome and the second being coordinate and the third being SNP names such as rs ID number. With this option, GTFtools will search for and output cis-SNPs for each gene annotated in the provided GTF file.<BR/>Example use:<BR/>gtftools -p snp_list.txt demo.gtf > cisSNP.bed||
| -f | -f specifies the upstream and downstream distance used to calculate cis-range of a gene. -f is specified in the format of 'distup-distdown', where distup represent the upstream distance from TSS and distdown means the downstream distance from the end of the gene. Note that this parameter takes effect only when the '-g' option is used. For example, using 'gtftools -g gene.bed -f 2000-1000 demo.gtf' means that 2000 bases upstream and 1000 bases downstream of the gene will be clculated as the cis-range and the cis-range will be output to the gene.bed file. By default, -f is set to 0-0, indicating that cis-range will not be calculated when using -g to calculate gene information.<BR/>Example use:<BR/>gtftools -g gene.bed -f 2000-1000 demo.gtf||
| -s | Output isoform coordination and parent-gene IDs in bed format.<BR/>Example use:<BR/>gtftools -s isoform.bed demo.gtf||
| -q | output 5' and 3' splice site regions in bed format. The region is based on MaxEntScan: the 5' donor site is 9 bases long with 3 bases in exon and 6 bases in intron, and the 3' acceptor site is 23 bases long with 20 bases in the intron and 3 bases in the exon.<BR/>Example use:<BR/>gtftools -q splice_regions.bed demo.gtf||
| -e | Output exons in bed format.<BR/>Example use:<BR/>gtftools -e exons.bed demo.gtf||
| -i | Output introns in bed format.<BR/>Example use:<BR/>gtftools -i introns.bed demo.gtf||
| -b | Output intergenic regions in bed format.<BR/>Example use:<BR/>gtftools -b intergenic_regions.bed demo.gtf||
| -k | Calculate introns (part of introns) that overlap with exons of other isoforms. The output is in bed format.<BR/>Example use:<BR/>gtftools -k introns_ovlp_exons.bed demo.gtf||
| -u | Output UTRs in bed format.<BR/>Example use:<BR/>gtftools -u utr.bed demo.gtf||
| -t | Output transcription start site (TSS)-flanking regions in bed format. It is calculated as (TSS-wup,TSS+wdown), where wup is a user-specified distance, say 1000bp, upstream of TSS, and wdown is the distance downstream of TSS. wup and wdown is defined by the w parameter specified by '-w'<BR/>Example use:<BR/>gtftools -t tss_regions.bed -w 1000-300 demo.gtf| Used for scanning TF-binding sites. |
| -w | w specifies the upstream and downstream distance from TSS as described in '-t'. w is specified in the format of wup-wdown, where wup and wdown represent the upstream and downstream distance of TSS. Default w = 1000-300 (that is, 1000 bases upstream of TSS and 300 bases downstream of TSS). This range is based on promoter regions used in the dbSNP database based on ref: Genome-wide promoter extraction and analysis in human, mouse, and rat, Genome Biology, 2005.<BR/>Example use:<BR/>gtftools -t tss_regions.bed -w 1000-300 demo.gtf||
| -c | Specify chromosomes to analyze. Dash(-) and comma(,) are allowed. For example, ‘-c 1-5,X,Y’ indicates 7 chromosomes: 1 to 5 together with X and Y. Default is 1-22, X and Y.<BR/>Example use (output genes which are on chromosomes 1, 2, X and Y):<BR/>gtftools -g gene.bed -c 1-2,X,Y demo.gtf<BR/>Example use (output splice site regions which are on chromosomes 1 and Y):<BR/>gtftools -q splice_regions.bed -c 1,Y demo.gtf||
| -v | Show version.<BR/>Example use:<BR/>gtftools -v||

## Help

In general, you can run 'gtftools -h' to obtain help documents.
Examples on how to use GTFtools are shown in the above Table.



## Contact

If any questions, please contact my tutor at:
Hongdong Li, hongdong@csu.edu.cn
