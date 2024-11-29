#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --account=ChrXh2
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=32G

cd /home/isdurso/ChrXh2/shannon/data/1000G
	
# Combine iPSYCH2015 autosomal data and 1000 Genomes for a PCA
# Convert 1000G pfile to binary format
# Exclude multi-allelic variants, and SNPs on MT, PAR1, PAR2, X and Y
# Only keep the variants present in iPSYCH2015
# Exclude strand ambiguous SNPs - if left in, these can cause a shift in the PC plots

cut -f1,2,3 /faststorage/jail/project/ChrXh2/data/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn.bim | sort > iPSYCH2015_snps.txt

plink2 --pfile all_phase3_ns \
	--max-alleles 2 \
	--autosome \
	--extract iPSYCH2015_snps.txt \
	--make-bed \
	--out 1000G

# Now attempt to merge the datasets
plink --bfile /faststorage/jail/project/ChrXh2/data/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn \
	--bmerge 1000G \
	--make-bed \
	--out iPSYCH2015_1000G_merge1
	
# This will produce multiple position warnings for several SNPs
grep -w Multiple iPSYCH2015_1000G_merge1.log | awk '{gsub(/['\''\.]/, "", $7); print $7}' > exclude_multi_pos_snps.txt

# Exclude variants with 3+ alleles
cat iPSYCH2015_1000G_merge1-merge.missnp exclude_multi_pos_snps.txt > exclude_multi_pos_allele.txt

plink2 --bfile 1000G \
	--exclude exclude_multi_pos_allele.txt \
	--make-bed \
	--out /scratch/$SLURM_JOB_ID/1000G_2

# LD prune the iPSYCH set before merging to speed up the computation

# Exclude the strand ambiguous SNPs (since I have not checked for strand alignment) and the long range ld regions
awk '(($5=="A" && $6=="T") || ($5=="T" && $6=="A") || ($5=="C" && $6=="G") || ($5=="G" && $6=="C"))' /faststorage/jail/project/ChrXh2/data/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn.bim  | awk '{print $2}' > ambig_snps.txt
cat ambig_snps.txt ../../scripts/long_range_ld_regions.txt > ambig_long_range_ld_regions.txt

plink --bfile /faststorage/jail/project/ChrXh2/data/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn \
	--snps-only \
	--maf 0.05 \
	--indep-pairwise 5000 300 0.05 \
	--exclude ambig_long_range_ld_regions.txt \
	--extract iPSYCH2015_snps.txt \
	--allow-no-sex \
	--out pruned_autosomes_iPSYCH2015_1000G_overlap
	
plink --bfile /faststorage/jail/project/ChrXh2/data/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn \
	--bmerge /scratch/$SLURM_JOB_ID/1000G_2 \
	--exclude exclude_multi_pos_allele.txt \
	--extract pruned_autosomes_iPSYCH2015_1000G_overlap.prune.in \
	--allow-no-sex \
	--make-bed \
