# CAGE_peak_annotation
A simple tool for CAGE peaks annotation

Makes annotation of CAGE peaks fast and simple.
Required input is BED file for CAGE peaks (-c) predicted by any available tool like DPI (https://github.com/hkawaji/dpi1), PromoterPipeline, CAGEr, etc. Column 7 must be a TSS position. Second required input is a file containing paths to UCSC tables of gene models (-u) and/or a file with paths to any custom BED files (-b). Final required option is an output folder (-o).

Other optional parameters are:
-t path to TEMP directory, where all temprary files will be saved;
-l length of desired extention regions (50 and 500 always used by default); several numbers could be provided, for example -l 200,1000,1500;
-a if Y/y the script will create input files for TssClassifier (https://sourceforge.net/p/tometools/wiki/TssClassifier/), default is N.

-u and -b files could be extended by additional column providing paths to tables with transcript and gene names (two columns):

<p>XM_025152649.1 LOC112532827</p>
<p>XM_025152753.1 LOC107049475</p>
<p>XM_025152751.1 LOC107049475</p>

Example folder includes all required and optional files for chicken CAGE peaks annotation and related output tables.

To make gene name table from UCSC annotation:
> cut -f2,13 ./example/ncbiRefSeq.txt > ./example/ncbiRefSeqToGeneName.txt

If custom BED table already includes gene names, please provide gene table anyway:
> awk 'BEGIN{OFS="\t"};{print $4, $4}'  ./example/Entrez_gene_galGal6.bed > ./example/Entrez_gene_galGal6_names.txt

Input -u and -b files:
> echo -e ./example/augustusGene.txt'\n'./example/ensGene.txt'\t'./example/ensemblToGeneName.txt'\n'./example/ncbiRefSeq.txt'\t'./example/ncbiRefSeqToGeneName.txt'\n'./example/genscan.txt > ./example/ucsc_path.txt
> echo -e ./example/Entrez_gene_galGal6.bed'\t'./example/Entrez_gene_galGal6_names.txt > ./example/bed_path.txt

If you use PromoterPipeline, BED file could be created by:
> sed '/^#/ d' level2.osc | tail -n +2 | awk 'BEGIN {OFS="\t"};{print $2, $3, $4, $1, $7,$5, $6}' > level2.bed

To run the script with example data try:
> ./CAGE_peaks_annotation.sh -c ./example/rDPI_merged_chicken6.bed -b ./example/bed_path.txt -u ./example/ucsc_path.txt -o ./example/output -t ./example/temp -l 1000 -a Y

Good luck! have fun :D
<p>Best,</p>
<p>Ruslan</p>

ruselusalbus@gmail.com | RMDevyatiyarov@kpfu.ru


