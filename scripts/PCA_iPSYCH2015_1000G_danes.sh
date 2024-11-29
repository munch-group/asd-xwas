#!/bin/bash
#SBATCH --job-name=PCA
#SBATCH --account=ChrXh2
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=32G

cd /home/isdurso/ChrXh2/shannon/data

# Create a file containing individuals to keep (1000G and iPSYCH)
# iPSYCH contains individuals whose grandparents were all born in Denmark
awk '{print $1, $2}' /faststorage/jail/project/ChrXh2/data/pca/ASD2015_pca2_ell_dk-dim3_8_8_8.menv.mds_cov > keep_danes.txt
awk '{print $1, $2}' /home/isdurso/ChrXh2/shannon/data/1000G/1000G.fam > keep_1000G.txt
cat keep_danes.txt keep_1000G.txt > keep_1000G_danes.txt

# Keep only overlapping SNPs
awk 'NR==FNR {seen[$2]; next} $2 in seen' 1000G/iPSYCH2015_1000G_pruned.bim 1000G/1000G.bim | awk '{print $2}' > 1000G/iPSYCH_1000G_pruned_overlapping_snps.txt

# Conduct PCA on combined data (iPSYCH2015 + 1000G)
# Use the plink2 --pca approx flag to conduct a PCA on > 5000 indiviuals (default = 10 PCs)
# Use the pruned autosomal data without ambiguous SNPs
plink2 --bfile 1000G/iPSYCH2015_1000G_pruned \
  --allow-no-sex \
  --keep keep_1000G_danes.txt \
  --extract 1000G/iPSYCH_1000G_pruned_overlapping_snps.txt \
  --pca approx \
  --out PCA_iPSYCH2015_1000G_danes
