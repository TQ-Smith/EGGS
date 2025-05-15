convertf -p param.txt
plink --file zeng2025 --recode vcf --out zeng2025
bgzip zeng2025.vcf
tabix zeng2025.vcf.gz
bcftools view -r 22 -e 'COUNT(GT="mis")=N_SAMPLES' zeng2025.vcf.gz | bcftools view -e 'N_ALT=0' | gzip > zeng2025_chr22.vcf.gz
