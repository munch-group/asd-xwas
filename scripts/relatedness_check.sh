#!/bin/bash
#SBATCH --job-name=ibd
#SBATCH --account=ChrXh2
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=16G

cd /home/isdurso/ChrXh2/shannon/data

# Exclusion list of SNPs in the MHC (chr 6, 25-35 MB) and CHR 8 inversion (chr 8, 7-13 MB)
awk '$1 == 6 && $4 >= 25000000 && $4 <= 35000000' /faststorage/jail/project/ChrXh2/data/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn.bim > /scratch/$SLURM_JOB_ID/tmp
awk '$1 == 8 && $4 >= 7000000 && $4 <= 13000000' /faststorage/jail/project/ChrXh2/data/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn.bim > /scratch/$SLURM_JOB_ID/tmp2
cat /scratch/$SLURM_JOB_ID/tmp /scratch/$SLURM_JOB_ID/tmp2 > long_range_ld_regions.txt

####
# Individual inclusion lists
# limit to those with asd phenotype available
awk '{print $1, $2}' /faststorage/jail/project/ChrXh2/data/pheno/asdGWAS2015.pheno > /scratch/$SLURM_JOB_ID/asd_phenotypes_available

# Subset to just those in the covariate file
# cov contains FID, IID, SOL?, C1 to C20, st1
awk '{print $1, $2}' /faststorage/jail/project/ChrXh2/data/pca/ASD2015_pca2_ell_dk-dim3_8_8_8.menv.mds_cov > /scratch/$SLURM_JOB_ID/IDs_asd_cov.txt

#####
cd /home/isdurso/ChrXh2/shannon/results

# Only keep those with phenotypes available
plink --bfile /faststorage/jail/project/ChrXh2/data/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn \
	--keep /scratch/$SLURM_JOB_ID/asd_phenotypes_available \
	--snps-only \
	--maf 0.05 \
	--indep-pairwise 5000 300 0.05 \
	--exclude ../scripts/long_range_ld_regions.txt \
	--allow-no-sex \
	--out ../data/pruned_autosomes
	
# Limit to individuals in the cov file (passing ancestry and relatedness checks)
plink --bfile /faststorage/jail/project/ChrXh2/data/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn \
	--extract ../data/pruned_autosomes.prune.in \
	--keep /scratch/$SLURM_JOB_ID/IDs_asd_cov.txt \
	--allow-no-sex \
	--make-bed \
	--out /scratch/$SLURM_JOB_ID/tmp

# Only report pairs with a minimum relatedness of 0.2 (limit output size)
plink --bfile /scratch/$SLURM_JOB_ID/tmp \
	--allow-no-sex \
	--genome \
	--min 0.2 \
	--out ibd_autosome
