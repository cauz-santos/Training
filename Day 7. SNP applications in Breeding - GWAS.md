## Week 2: Applied Bioinformatics for Genomics and Breeding

## Day 7: SNP Applications in Breeding - GWAS and Gene Interpretation

Welcome to Day 7!

Today we connect **SNP genotypes** to a **quantitative trait** (sucrose content, **SUC**) using **Genome-Wide Association Studies (GWAS)**, and then interpret significant signals biologically by mapping top SNPs to **genes** and **functions**.

This is a **hands-on** session; you will:
- Prepare genotype/phenotype/covariates for GWAS
- Run GWAS in **PLINK** (with population-structure correction via **PC covariates** from Day 6)
- Visualize results (**Manhattan** and **QQ** plots in R)
- Identify **top SNPs**, map them to **genes** (GFF3), and extract **putative functions**
- Discuss how hits feed into **Marker-Assisted Selection (MAS)** and **Genomic Selection (GS)**

> **Relevance:**  
> GWAS uncovers **genomic regions** that control traits like sucrose. These become **markers** for MAS, and inform **genome-wide prediction** in GS—shortening cycles and improving accuracy in breeding pipelines.

### Input Data

- `dataset120_chr18.vcf.gz` — high-quality, biallelic SNPs
- `pca_results.eigenvec` — PCA eigenvectors (Day 6 output)  
- `gwas_phen_table_120.csv` — phenotypes with **SUC** column (g/100g dry matter)

**Example `phenotypes.csv` (CSV with header):**
```
IID,SUC
PDAC253_Ram_S,17.88
PDAC254_Rij_S,20.86
PDAC255_Sin_S,21.24
...
```

> **Assumptions:**  
> - **IID** matches the **second column** of your PLINK `.fam` file (IID).  
> - You can set **FID = IID** if family IDs are not used.  
> - Missing phenotypes are `NA`.

---

### Part 0 — Prepare Genotypes, Phenotypes, and Covariates

First create a folder for the file outputs of day 7 in your home directory:
   ```bash
   mkdir 07_gwas_selection
   ```
Now, enter the folder with the command `cd`

Then create some subfolders, for each specific analysis we will perform:
   ```bash
   mkdir plink gwas mas gs-lite
   ```


### Step 1 — Convert VCF → PLINK (binary)

Create and submit:

```bash
vi 00_vcf2plink_gwas.sh
```
copy and paste:
```bash
#!/bin/bash
#SBATCH --job-name=vcf2plink_gwas
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH -o vcf2plink_gwas.out
#SBATCH -e vcf2plink_gwas.err

module load PLINK

IN_VCF="/lisc/scratch/course/pgbiow/data/VCF/dataset120_chr18.vcf.gz"

echo "Converting VCF to PLINK binary..."
plink --vcf "${IN_VCF}" \
      --make-bed \
      --allow-extra-chr \
      --double-id \
      --out ./plink/gwas_data

echo "Done: gwas_data.bed/.bim/.fam"
```
save the file and then submit the job:
```bash
sbatch 00_vcf2plink_gwas.sh
```

**Outputs:** `gwas_data.bed/.bim/.fam`


### Step 2 — Build a PLINK-friendly phenotype file for **SUC**

We’ll construct **`pheno_suc.txt`** with 3 columns: `FID IID SUC`.

```bash
vi 01_make_pheno.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=make_pheno
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH -o make_pheno.out
#SBATCH -e make_pheno.err

set -euo pipefail

PHENO_CSV="/lisc/scratch/course/pgbiow/data/metadata/gwas_phen_table_120.csv"   # columns: IID,SUC
FAM="./plink/gwas_data.fam"

# Create header
echo -e "FID\tIID\tSUC" > pheno_suc.txt

# Build associative array from phenotypes (CSV -> TSV)
awk -F, 'NR>1{ph[$1]=$2} END{for(k in ph) print k"\t"ph[k] }' "${PHENO_CSV}" \
| sort -k1,1 > pheno_suc_tmp.tsv

# Append FID IID SUC (FID=IID if FID unknown)
awk 'NR==FNR{ph[$1]=$2; next} {fid=$1; iid=$2; val=(iid in ph?ph[iid]:"NA"); print fid"\t"iid"\t"val}' \
    pheno_suc_tmp.tsv "${FAM}" >> pheno_suc.txt

rm -f pheno_suc_tmp.tsv

echo "Done: pheno_suc.txt (FID IID SUC)"
```

```bash
sbatch 01_make_pheno.sh
```

### Step 3 — Prepare **PCA covariates** (PC1–PC5) from Day 6
Add PC1–PC5 as covariates so GWAS controls for broad genetic background (ancestry/relatedness) differences among samples. Without them, some SNPs look “significant” just because certain groups both share those alleles and tend to have different sucrose values. PCs capture that group effect; adjusting for them removes this bias and leaves signals that are more likely to be truly linked to sucrose.

Create **`covar_pcs.txt`** with header `FID IID PC1 PC2 PC3 PC4 PC5`:

```bash
vi 02_make_covariates.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=make_covariates
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH -o make_covariates.out
#SBATCH -e make_covariates.err

EIGENVEC="../06_diversity_structure/pca/pca_results.eigenvec"

# pca_results.eigenvec: FID IID PC1 PC2 ...
# Keep first 5 PCs (adjust if needed)
{ 
  echo -e "FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5"
  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' "${EIGENVEC}"
} > covar_pcs.txt

echo "Done: covar_pcs.txt (PC1..PC5)"
```

```bash
sbatch 02_make_covariates.sh
```

> **Why PCs as covariates?**  
> PCs capture **population structure**; including them reduces **false positives** (spurious associations due to stratification).

---

### Part 1 — GWAS QC (basic genotype-level filters)

We’ll apply minimal QC commonly used before association:

- **MAF filter**: remove very rare variants (e.g., MAF < 0.05)  
- **Missingness per SNP**: remove poorly genotyped sites (e.g., >20% missing)  

> We apply MAF, missingness, and HWE filters to remove rare, error-prone, or inconsistent SNPs so GWAS tests high-quality markers, improving power and reducing false positives.

```bash
vi 10_qc_make_subset.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=gwas_qc
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:40:00
#SBATCH -o gwas_qc.out
#SBATCH -e gwas_qc.err

module load PLINK

plink --bfile ./plink/gwas_data \
      --maf 0.05 \
      --geno 0.20 \
      --make-bed \
      --allow-extra-chr \
      --out ./plink/gwas_data_qc

echo "QC subset: gwas_data_qc.*"
```

```bash
sbatch 10_qc_make_subset.sh
```

> **Note:** Thresholds are dataset-dependent. For training we use pragmatic defaults.

---

### Part 2 — GWAS for **SUC** with PLINK (linear model + PC covariates)

We’ll run an **additive linear regression**, including **PC1–PC5** as covariates.

**What is “additive linear regression” (for our plants)?**

In GWAS, **additive linear regression** tests whether the trait changes in a straight-line way with the **number of alternate alleles** at a SNP for each **individual/plant**:
- Genotype is coded **0, 1, or 2** (copies of the alternate allele).
- We fit: `Trait ≈ Intercept + BETA × Genotype (+ covariates like PC1–PC5)`.
- **BETA** is the **per-allele effect**: how much the trait (e.g., sucrose, g/100g) changes **for each extra copy** of the alternate allele.

**Example:** If `BETA = +0.8`, then (on average)
- Genotype **0** → baseline SUC  
- Genotype **1** → baseline **+ 0.8**  
- Genotype **2** → baseline **+ 1.6**

This assumes **additivity** (heterozygotes are roughly halfway between the two homozygotes). If you suspect dominance/overdominance, you’d use a genotypic or dominance model—but the additive model is the standard, efficient default for genome-wide scans.

```bash
vi 20_run_gwas_suc.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=gwas_suc
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH -o gwas_suc.out
#SBATCH -e gwas_suc.err

module load PLINK

plink --bfile ./plink/gwas_data_qc \
      --pheno pheno_suc.txt \
      --pheno-name SUC \
      --covar covar_pcs.txt \
      --covar-name PC1,PC2,PC3,PC4,PC5 \
      --allow-extra-chr \
      --linear hide-covar \
      --allow-no-sex \
      --out ./gwas/gwas_suc_linear

echo "Done: gwas_suc_linear.assoc.linear"
```

```bash
sbatch 20_run_gwas_suc.sh
```

**Output:** `gwas_suc_linear.assoc.linear`  
- Contains columns: `CHR BP SNP A1 TEST NMISS BETA STAT P ...`  
- Use rows with `TEST == "ADD"` for additive model results.

> **Tip:** If your trait is **non-normal**, consider rank-normalizing phenotypes or using robust models. For binary traits use `--logistic`.
>
> **Relevance of GWAS**
> - Identifies **trait-linked markers** that can be turned into **KASP/SNP-chip assays**.
> - Guides **parent selection** and **cross design** by highlighting favorable alleles/haplotypes.
> - Prioritizes **validation targets** to de-risk downstream R&D spending.
> - Feeds **genomic selection** models with signal-rich loci, improving prediction accuracy.
---

### Part 3 — Visualize GWAS (Manhattan + QQ) in R

Now we will load the PLINK GWAS results and create two standard checks:  
- **Manhattan plot:** each point is a SNP; the y-axis is **−log10(p)**, so taller points mean stronger associations. Peaks across the x-axis (chromosomes) highlight regions linked to **SUC**.  
- **QQ plot:** compares observed p-values to what we’d expect by chance; big upward deviations suggest inflation or true signals.  
We’ll also add a **Bonferroni** genome-wide line (cutoff = 0.05 / number of SNPs) to control false positives, and compute **FDR (BH)** as a less conservative alternative.

Open an **interactive R** session (with graphics) on the cluster or RStudio.

```r
# ================================
# GWAS Visualization: SUC (Manhattan plot)
# ================================

library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- Load PLINK linear results ---
dt <- fread("./gwas/gwas_suc_linear.assoc.linear")
table(dt$TEST, useNA = "ifany")

# --- Keep only additive model and valid P-values ---
dt <- dt[TEST == "ADD"]
dt[, P := as.numeric(P)]
dt <- dt[!is.na(P)]

# --- Fix CHR column: extract numeric index from labels like "chrLG1", "chrLG2", etc. ---
dt[, CHR := as.numeric(gsub("chrLG", "", CHR))]
dt <- dt[!is.na(CHR) & !is.na(BP)]

# --- Genomic inflation (lambda) ---
dt[, CHISQ := qchisq(1 - P, df = 1)]
lambda <- median(dt$CHISQ, na.rm = TRUE) / 0.456
cat("Genomic inflation factor (lambda):", round(lambda, 3), "\n")

# --- Multiple testing thresholds ---
M <- nrow(dt)
bonf <- 0.05 / M
dt[, FDR_BH := p.adjust(P, method = "BH")]
cat("Bonferroni threshold (alpha = 0.05):", signif(bonf, 3), "\n")

# --- Cumulative position for Manhattan plot ---
setorder(dt, CHR, BP)
chr_sizes <- dt[, .(chr_len = max(BP)), by = CHR]
chr_sizes[, chr_start := cumsum(shift(chr_len, fill = 0))]
dt <- merge(dt, chr_sizes[, .(CHR, chr_start)], by = "CHR", all.x = TRUE)
dt[, BPcum := BP + chr_start]

# --- Axis positions (center of each chromosome block) ---
axis_df <- dt[, .(center = (min(BPcum) + max(BPcum)) / 2), by = CHR]

# --- Manhattan plot ---
dt[, logP := -log10(P)]
dt[, color := as.factor(CHR %% 2)]

p_manhattan <- ggplot(dt, aes(x = BPcum, y = logP, color = color)) +
  geom_point(alpha = 0.75, size = 1) +
  scale_color_manual(values = c("steelblue", "darkgrey")) +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", color = "red", linewidth = 0.7) +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Chromosome", y = "-log10(p-value)",
       title = "GWAS Manhattan Plot — Sucrose (SUC)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Preview plot
print(p_manhattan)

# --- QQ plot ---
observed <- sort(dt$P)
expected <- -log10(ppoints(length(observed)))
observed_logp <- -log10(observed)

qq_df <- data.frame(Expected = expected, Observed = observed_logp)

p_qq <- ggplot(qq_df, aes(x = Expected, y = Observed)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = sprintf("GWAS QQ Plot (SUC) — lambda = %.3f", lambda),
       x = "Expected -log10(p)",
       y = "Observed -log10(p)") +
  theme_minimal(base_size = 14)

# Preview plot
print(p_qq)

# --- Save plots to PDF after inspection ---
ggsave("./gwas/GWAS_SUC_Manhattan.pdf", plot = p_manhattan, width = 10, height = 5)
ggsave("./gwas/GWAS_SUC_QQ.pdf", plot = p_qq, width = 6, height = 6)

# --- Export top hits ---
setorder(dt, P)
fwrite(dt[1:20], "./gwas/top20_hits_SUC.tsv", sep = "\t")
fwrite(dt[P <= bonf], "./gwas/bonferroni_hits_SUC.tsv", sep = "\t")

cat("PDFs saved: GWAS_SUC_Manhattan.pdf, GWAS_SUC_QQ.pdf\n")
cat("Tables saved: top20_hits_SUC.tsv, bonferroni_hits_SUC.tsv\n")
```

### Inspecting GWAS Bonferroni Hits (SUC)

After running the analysis, a TSV file with all SNPs passing the Bonferroni correction (α = 0.05) is saved at:

```bash
cd ~/07_gwas_selection/gwas

# View the first 10 lines
head bonferroni_hits_SUC.tsv

# Scroll interactively
less bonferroni_hits_SUC.tsv

# View just SNP and p-value columns
cut -f2,9 bonferroni_hits_SUC.tsv | head
```

This file contains a list of SNPs that are **significantly associated with sucrose (SUC) concentration** based on the additive model (`TEST == "ADD"`). Each row represents one SNP that passed the Bonferroni threshold for genome-wide significance.

#### What’s in the file?

| Column         | Description |
|----------------|-------------|
| `CHR`          | Chromosome number (parsed from the original `chrLG` format) |
| `SNP`          | SNP identifier, usually in the format `chr:position:alleles` |
| `BP`           | Base-pair position on the chromosome |
| `A1`           | Effect allele |
| `TEST`         | Model used in GWAS (here, always `ADD`) |
| `NMISS`        | Number of non-missing genotypes used in the association test |
| `BETA`         | Estimated effect size of the allele on the phenotype |
| `STAT`         | Test statistic (e.g., t-score) |
| `P`            | Raw p-value for the association |
| `CHISQ`        | Chi-square value computed from the p-value |
| `FDR_BH`       | Adjusted p-value using Benjamini-Hochberg False Discovery Rate |
| `chr_start`    | Cumulative base-pair offset used for Manhattan plotting |
| `BPcum`        | Cumulative base-pair position for plotting across chromosomes |
| `logP`         | -log10(P), used for plotting |
| `color`        | Used internally to alternate Manhattan plot colors (e.g., 0/1) |


**Interpretation checks:**  
- Are there **clear peaks** in the Manhattan plot?  
- Is the **QQ plot** close to expectation (no heavy inflation)?  
- Do top hits remain after **FDR** or **Bonferroni** correction?

> **Relevance:**  
> Peaks identify **candidate loci** controlling sucrose—targets for **marker development**, **introgression**, and **GS validation**.


### Regional Association Plot (±250 kb) around the Top Hit

Now we’ll zoom into the strongest GWAS peak to inspect the local signal shape (is it a tight single peak or a broad region?).

**Do I need to change `chr` and `pos`?**  
**No** for the automatic version below — it **reads the top SNP from your PLINK results** and uses its `CHR`/`BP`.  
Only change them if you want to zoom into a **specific SNP or region** (manual variants provided).

```r
# --- Regional association plot around the top SNP (±250 kb) ---

library(data.table)
library(ggplot2)
library(ggrepel)

# Load additive model GWAS results
dt <- fread("gwas/gwas_suc_linear.assoc.linear")[TEST == "ADD" & !is.na(P)]

# Extract numeric chromosome from "chrLG14" format
dt[, CHR_NUM := as.numeric(gsub("chrLG", "", CHR))]

# Identify top SNP (lowest p-value)
setorder(dt, P)
lead <- dt[1]
chr  <- lead$CHR_NUM
pos  <- lead$BP

cat(sprintf("Top SNP: %s (CHR %s, POS %d, P = %.2e)\n", lead$SNP, chr, pos, lead$P))

# Subset ±250 kb window around top SNP
win <- dt[CHR_NUM == chr & BP >= (pos - 250000) & BP <= (pos + 250000)]

# Identify top SNP in the window
top_win <- win[which.min(P)]

# --- Plot to screen ---
p_regional <- ggplot(win, aes(x = BP, y = -log10(P))) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = pos, linetype = "dashed", color = "red") +
  geom_text_repel(
    data = top_win,
    aes(label = SNP),
    size = 3,
    nudge_y = 0.5,
    box.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  labs(
    title = sprintf("Regional Association: Chr %s ±250 kb around %s", chr, format(pos, big.mark = ",")),
    x = "Genomic Position (bp)",
    y = "-log10(P)"
  ) +
  theme_minimal(base_size = 14)

# Show the plot
print(p_regional)

# --- Save to file after inspection ---
ggsave("GWAS_SUC_RegionalTop_labeled.png", plot = p_regional, width = 10, height = 4, dpi = 150)

cat("✅ Saved: GWAS_SUC_RegionalTop_labeled.png\n")
```
---

### Part 4 — From SNP to Gene (annotation & function)

Goal: For our **top SNP(s)**, find the **overlapping/nearest gene(s)** and a **putative function**.

> **You need:**  
> - The **same reference genome** used for alignment & variant calling (FASTA)  
> - Its **annotation** in **GFF3** (gene coordinates + attributes, e.g., `ID`, `Name`, `product`)  
> - **BEDTools** module (for genomic intersections)

---

### 📂 Files for Date Palm Genome

We are using the **Date Palm genome** from NCBI  
👉 [GCA_009389715.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009389715.1/)

These files are available on the cluster:  
`/lisc/scratch/course/pgbiow/data/genomes/`  
date_palm_genome.fna # reference genome (FASTA)  
date_palm_genomic.gff # annotation in GFF3 format  
date_palm_genomic.gtf # annotation in GTF format (alternative)  
 
**Make a BED of top SNPs:**  

We’ll map **Bonferroni-significant** hits (or top 20 if none pass) to genes. This script creates:  
- `top_snps.bed` — BED coordinates of your strongest GWAS SNPs (Bonferroni-significant or top 20).  
- `genes.bed` — BED coordinates of all genes parsed from your GFF3 with IDs/names/products.  
You’ll use these with BEDTools to find which genes overlap or are nearest to the top SNPs.

```bash
vi 30_prepare_bed_for_annotation.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=prep_annot
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH -o prep_annot.out
#SBATCH -e prep_annot.err

set -euo pipefail

# --- Reference annotation (Date Palm genome from NCBI) ---
GENOME_GFF="/lisc/scratch/course/pgbiow/data/genomes/date_palm_genomic.gff"

# --- Create output folder ---
OUTDIR="annotation"
mkdir -p $OUTDIR

# --- Step 1) Pick SNPs to annotate (Bonferroni first, else top20) ---
if [ -s bonferroni_hits_SUC.tsv ]; then
  echo "Using Bonferroni-significant SNPs..."
  awk 'NR>1 {print $1, $3, $2, $9}' OFS='\t' bonferroni_hits_SUC.tsv > $OUTDIR/snps_for_annot.tsv
else
  echo "No Bonferroni file, using Top20 hits..."
  awk 'NR>1 {print $1, $3, $2, $9}' OFS='\t' top20_hits_SUC.tsv > $OUTDIR/snps_for_annot.tsv
fi
# Output columns: CHR BP SNP P

echo "Wrote: $OUTDIR/snps_for_annot.tsv"

# --- Step 2) Convert SNPs into 0-based BED format ---
awk 'BEGIN{OFS="\t"} {start=$2-1; end=$2; print $1, start, end, $3}' $OUTDIR/snps_for_annot.tsv > $OUTDIR/top_snps.bed

echo "Wrote: $OUTDIR/top_snps.bed"

# --- Step 3) Extract gene features from GFF3 into BED ---
awk -v OFS="\t" '$3=="gene" {
  split($9,a,";"); id="NA"; name="NA"; prod="NA";
  for(i in a){
    if(a[i] ~ /^ID=/){sub(/^ID=/,"",a[i]); id=a[i]}
    if(a[i] ~ /^Name=/){sub(/^Name=/,"",a[i]); name=a[i]}
    if(a[i] ~ /^product=/){sub(/^product=/,"",a[i]); prod=a[i]}
    if(a[i] ~ /^gene=/){sub(/^gene=/,"",a[i]); name=a[i]} # fallback if only "gene=" exists
  }
  print $1, $4-1, $5, id "|" name "|" prod
}' "$GENOME_GFF" > $OUTDIR/genes.bed

echo "Wrote: $OUTDIR/genes.bed (gene features)"

# --- Step 4) Intersect SNPs with genes ---
module load BEDTools/2.30.0  # adjust module name if different on your cluster

bedtools intersect -a $OUTDIR/top_snps.bed -b $OUTDIR/genes.bed -wa -wb > $OUTDIR/snps_in_genes.tsv
bedtools closest   -a $OUTDIR/top_snps.bed -b $OUTDIR/genes.bed > $OUTDIR/snps_nearest_genes.tsv

echo "Wrote: $OUTDIR/snps_in_genes.tsv (overlaps)"
echo "Wrote: $OUTDIR/snps_nearest_genes.tsv (closest genes)"
```

```bash
sbatch 30_prepare_bed_for_annotation.sh
```

**Intersect SNPs with genes (overlap) and find nearest genes:**

We’ll use **BEDTools** for exact overlaps and the **nearest** gene within ±10 kb (change as needed).

```bash
vi 31_map_snps_to_genes.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=map_snps_genes
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH -o map_snps_genes.out
#SBATCH -e map_snps_genes.err

module load bedtools

# Overlapping genes
bedtools intersect -a top_snps.bed -b genes.bed -wa -wb > snp_gene_overlaps.tsv

# Nearest gene (report distance)
bedtools closest -a top_snps.bed -b genes.bed -d > snp_gene_nearest.tsv

echo "Created: snp_gene_overlaps.tsv (overlaps), snp_gene_nearest.tsv (nearest with distance)"
```

```bash
sbatch 31_map_snps_to_genes.sh
```

**Interpreting outputs:**
- `snp_gene_overlaps.tsv` → if a SNP lies **within** a gene (exon/intron span)  
- `snp_gene_nearest.tsv` → the **closest** gene and the **distance** (0 if overlapping)

**Summarize mappings (friendly table)**

We’ll produce an easy-to-read table linking **SNP → gene → putative function**.

```bash
vi 32_summarize_annotations.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=sum_annot
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH -o sum_annot.out
#SBATCH -e sum_annot.err

# If overlaps exist, prefer them; otherwise use nearest within 10kb
# top_snps.bed: chr start end snp
# genes.bed: chr start end id|name|product
# snp_gene_overlaps.tsv / snp_gene_nearest.tsv: a(4 cols) + b(4 cols) + [distance]

echo -e "SNP\tCHR\tBP\tGene_ID\tGene_Name\tProduct\tRelation\tDistance_bp" > snp_gene_summary.tsv

# 1) Overlaps
if [ -s snp_gene_overlaps.tsv ]; then
  awk 'BEGIN{OFS="\t"}{
    split($8,a,"|"); gid=a[1]; gname=a[2]; prod=a[3];
    bp=$2+1;
    print $4,$1,bp,gid,gname,prod,"overlap",0
  }' snp_gene_overlaps.tsv >> snp_gene_summary.tsv
fi

# 2) Nearest within 10000 bp (10 kb), excluding those already printed
# Build a set of SNPs already added
cut -f1 snp_gene_summary.tsv | tail -n +2 | sort -u > already.tsv

awk -vMAXD=10000 'BEGIN{OFS="\t"}
  FNR==NR {seen[$1]=1; next}
  {
    split($8,a,"|"); gid=a[1]; gname=a[2]; prod=a[3]; dist=$9
    snp=$4; chr=$1; bp=$2+1
    if(!(snp in seen) && dist<=MAXD){
      print snp, chr, bp, gid, gname, prod, "nearest", dist
      seen[snp]=1
    }
  }' already.tsv snp_gene_nearest.tsv >> snp_gene_summary.tsv

rm -f already.tsv

echo "Wrote: snp_gene_summary.tsv"
```

```bash
sbatch 32_summarize_annotations.sh
```

Open `snp_gene_summary.tsv` — you’ll see each top SNP, the overlapping/nearest gene, and a **product/description** (when present in GFF3 attributes).

> **If function is missing:**  
> - Your GFF annotation may not include `product`—try `Name`/`gene` attributes or consult the genome’s annotation README.  
> - For deeper function (domains/GO), run tools like **InterProScan** offline on the protein FASTA of candidate genes (advanced, not covered here).

> **Relevance:**  
> - Translating SNPs to **genes and functions** turns statistical signals into **biological hypotheses**, which alleles/genes to track, validate, and deploy in breeding.
> - Accelerates **marker development** and **IP differentiation** around candidate genes/alleles.
> - Informs **assay design** (haplotype tagging, primers) and prioritizes **functional validation** (expression, knockouts).
> - Creates a pipeline from **statistical hit → deployable marker → breeding action**.
---

### Part 5 — From GWAS to Selection Decisions

Now that we have GWAS hits and candidate genes, we’ll take the final step: show how these results translate into breeding tools through MAS (diagnostic markers) and GS-lite (genomic scores).

You will:  
1) build a **MAS marker table** from your top SNPs  
2) compute a quick **polygenic (GS-lite) score** to rank lines  
*Everything below assumes you already ran Part 3 and have: `gwas_suc_linear.assoc.linear`, `bonferroni_hits_SUC.tsv` and/or `top20_hits_SUC.tsv`, plus your genotype set `gwas_data_qc.*` and `phenotypes.csv` (IID,SUC).


**A) Make a MAS-ready marker table (effect allele, effect size, p-value)**  
We will build a script that collects the top SNPs from GWAS results (Bonferroni-significant or top 20), joins them with PLINK’s BIM file to add chromosome, position, and alleles, and outputs a ready-to-use marker list for marker-assisted selection (MAS).

```bash
vi 40_build_mas_markers.sh
```

Paste the following content:

```bash
#!/bin/bash
#SBATCH --job-name=mas_markers
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH -o mas_markers.out
#SBATCH -e mas_markers.err

set -euo pipefail

# 1) Pick the source of top hits (prefer Bonferroni, else Top20)
if [ -s bonferroni_hits_SUC.tsv ]; then SRC=bonferroni_hits_SUC.tsv; else SRC=top20_hits_SUC.tsv; fi

# 2) From PLINK BIM, build a lookup: SNP -> CHR BP A1 A2
#                 SNP  CHR  BP   A1   A2
awk 'BEGIN{OFS="\t"}{print $2, $1, $4, $5, $6}' gwas_data_qc.bim > bim.lookup

# 3) From the hits table, extract SNP, P, BETA by header names (keeps header first)
awk 'BEGIN{FS=OFS="\t"}
NR==1 {for(i=1;i<=NF;i++) h[$i]=i; print "SNP","P","BETA"; next}
      {print $h["SNP"], $h["P"], $h["BETA"]}' "$SRC" > assoc.slim

# 4) Join hits with BIM to add CHR/BP/A1/A2, then reorder columns and add a header
sort -k1,1 assoc.slim > a
sort -k1,1 bim.lookup > b
join -t $'\t' -1 1 -2 1 a b | awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$6,$7,$3,$2}' > mas_markers.tsv
(echo -e "SNP\tCHR\tBP\tA1(effect)\tA2\tBETA\tP"; cat mas_markers.tsv) > tmp && mv tmp mas_markers.tsv

# 5) (Optional) Also write a plain SNP list (no header)
cut -f1 mas_markers.tsv | tail -n +2 > mas_markers.snplist

# Clean temp
rm -f a b bim.lookup assoc.slim

echo "Created: mas_markers.tsv  (and mas_markers.snplist)"
```

Save and submit:

```bash
sbatch 40_build_mas_markers.sh
```

> Relevance: A clean list of deployable markers (effect allele, effect size, position) to design assays and screen parents/progeny.

### B) Quick GS-lite: Compute a Polygenic Score (PGS) and Rank Lines
After GWAS, we know which SNPs are significantly associated with sucrose and their estimated effect sizes (**BETA**). In **Marker-Assisted Selection (MAS)** we could take a few diagnostic SNPs. But in **Genomic Selection (GS)** we want to use many SNPs together to rank lines.  

Here we build a **Polygenic Score (PGS)**: for each plant, we multiply the genotype (0/1/2 copies of the allele) by the SNP’s effect size (BETA) and sum across top SNPs. This gives one number per line = the “genomic score.” If PGS correlates with observed sucrose, we can already prioritize top candidates. "Plants with higher scores should, on average, have higher sucrose. This is a mini version of Genomic Selection (GS)."

```bash
# 1) Build weights file for PLINK --score (columns: SNP  A1  BETA)
awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1,$4,$6}' mas_markers.tsv > gs_weights.tsv
echo "Wrote gs_weights.tsv (SNP A1 BETA)"

# 2) Compute per-line polygenic score (PGS) using PLINK
module load plink
plink --bfile gwas_data_qc \
      --score gs_weights.tsv 1 2 3 header sum \
      --out gs_pgs
# Output: gs_pgs.profile  (contains FID IID SCORE)

# 3) Rank lines by PGS and compare to observed SUC phenotype
awk 'NR==1{for(i=1;i<=NF;i++)if($i=="SCORE")s=i; next}{print $1"\t"$2"\t"$s}' gs_pgs.profile > pgs.tsv
awk -F, 'NR>1{print $1"\t"$2}' phenotypes.csv | sort -k1,1 > pheno.tsv
sort -k2,2 pgs.tsv > pgs.sorted
join -t $'\t' -1 2 -2 1 pgs.sorted pheno.tsv | awk 'BEGIN{OFS="\t"}{print $2,$1,$3,$4}' > pgs_with_suc.tsv
(echo -e "FID\tIID\tPGS\tSUC"; cat pgs_with_suc.tsv) > tmp && mv tmp pgs_with_suc.tsv
sort -k3,3gr pgs_with_suc.tsv > pgs_ranked.tsv
head -n 21 pgs_ranked.tsv > top20_by_PGS.tsv
echo "Created: pgs_ranked.tsv (all ranked) and top20_by_PGS.tsv"

# 4) Quick R plot: check correlation between PGS and SUC
R --vanilla <<'EOF'
d <- read.table("pgs_ranked.tsv", header=TRUE, sep="\t")
png("PGS_vs_SUC.png", 1000, 800, res=150)
plot(d$PGS, d$SUC, xlab="Polygenic Score (PGS)", ylab="Sucrose (SUC)", main="PGS vs SUC (quick check)")
abline(lm(SUC ~ PGS, data=d))
legend("topleft", bty="n", legend=paste0("r = ", round(cor(d$PGS, d$SUC, use='complete.obs'),3)))
dev.off()
EOF
echo "Saved: PGS_vs_SUC.png"
```

**Interpretation:**  
`pgs_ranked.tsv` → all lines ranked by genomic score.  
`PGS_vs_SUC.png` → scatterplot showing how well the genomic score predicts sucrose.  
If correlation (`r`) is positive and strong, even this “lite” version shows promise.  

> Relevance GS-lite:
> This gives a fast, interpretable genomic ranking you can already use today for pre-selection, while waiting for larger training sets and more advanced GS models. It shows how GWAS results can directly feed into breeding decisions.

## You have completed **Day 7**!

---

### Useful Tutorials and Resources

- [PLINK Association Analysis (GWAS)](https://zzz.bwh.harvard.edu/plink/anal.shtml)  
- [qqman R package (Manhattan & QQ plots)](https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html)
- [How to perform GWAS](https://frontlinegenomics.com/how-to-perform-a-genome-wide-association-study-gwas/)
- [GEMMA Documentation (Mixed Model GWAS)](https://github.com/genetics-statistics/GEMMA)  
- [GCTA fastGWA Manual](http://cnsgenomics.com/software/gcta/#GWAS)  

https://frontlinegenomics.com/how-to-perform-a-genome-wide-association-study-gwas/

### Appendix (Optional/Advanced): Mixed Models & When to Use Other GWAS Tools

For related/structured samples, prefer **mixed linear models (MLMs)** that include a **kinship matrix** to control sample relatedness and residual structure.
*Not required for this 4-hour practical, but recommended for production analyses.*

**When to use alternatives**
- **Strong relatedness / stratification:** use **GEMMA**, **GCTA fastGWA**, **EMMAX**, **BOLT-LMM**.  
- **Very large cohorts (speed):** **BOLT-LMM**, **GCTA fastGWA**.  
- **Imbalanced case–control / rare binary traits:** **SAIGE**.  
- **Plant-focused methods (MLM, FarmCPU, BLINK):** **GAPIT** (R), **TASSEL**.  
- **Fine-mapping after peaks:** **FINEMAP**, **SuSiE**, **CAVIAR**, **coloc**.
