## Week 2: Applied Bioinformatics for Genomics and Breeding

## Day 6: Diversity and Structure Analysis with SNPs

Welcome to Day 6 of our bioinformatics training!  

Today, we take a major step forward — from simply **generating variant data** to actually **understanding the biology behind the variation**. Using SNPs (Single Nucleotide Polymorphisms), we will explore how genomes differ between individuals and how those differences reflect patterns of **genetic diversity**, **population structure**, and **ancestry**.

This is a **hands-on session** designed to run on an HPC cluster using **Slurm**.


### Learning Objectives

By the end of this session, you will be able to:

- Understand the biological meaning behind genetic structure and diversity
- Prepare and convert VCF files for population genetics tools
- Perform and interpret **Principal Component Analysis (PCA)**
- Run **Admixture analysis** to infer individual ancestries
- Generate **SNP density plots** across chromosomes
- Estimate **basic population diversity statistics** (e.g., observed heterozygosity, missingness, allele frequencies)
- Use **Slurm job scripts** to automate analyses on a cluster
- Begin to critically evaluate and interpret population-level genomic data  



### Input Data
For training purposes, we will use a SNP dataset from Verdant:  

- **File**: `Report_DOp25-10208_4.1.vcf`  
- **Individuals**: **282 individuals**.
- **Variants**: **6878 SNPs**.  

This VCF enables us to run the full pipeline — including **PCA, ADMIXTURE, LD pruning/decay, and diversity statistics** — in a practical timeframe.


### Why Does This Matter?

Understanding **genetic diversity and structure** is essential in genomics and breeding:

- In **breeding programs**, it helps identify divergent individuals, avoid inbreeding, and guide parent selection.
- In **population genetics**, it reveals how populations are related and how they evolved.
- In **conservation**, it supports decisions on which populations or individuals to prioritize.

The tools we explore today — like PCA and ADMIXTURE — are **foundational methods** in modern genomics. Combined with diversity metrics (e.g., heterozygosity, allele frequency spectra), they allow you to **quantify and visualize genetic patterns** in your dataset.


**Table 1. Tools You Will Use**

| Tool       | Purpose                                      |
|------------|----------------------------------------------|
| **PLINK**      | Data conversion, PCA, diversity stats        |
| **ADMIXTURE**  | Ancestry inference (admixture proportions)  |
| **VCFtools**   | Summary statistics, missingness, MAF        |
| **BCFtools**   | SNP extraction, VCF querying                 |
| **R / Python** | Visualization (PCA plots, admixture bars, SNP density) |

---

## Part 0 — Data Preparation (VCF → PLINK + LD Pruning)

First create a folder for the file outputs of day 6 in your home directory:
   ```bash
   mkdir 06_diversity_structure
   ```
Now, enter the folder with the command `cd`

Then create some subfolders, for each specific analysis we will perform:
   ```bash
   mkdir plink diversity pca admixture snv_density
   ```

We will convert a **filtered VCF** into PLINK format and perform **LD pruning**. The pruned dataset will be used for **PCA** (and later for ADMIXTURE/GWAS).

### Step 1 — Clean the VCF (dedup + filters)  

Create a file called `00a_clean_vcf.sh`:

```bash
vi 00a_clean_vcf.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=clean_vcf
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH -o clean_vcf.out
#SBATCH -e clean_vcf.err

set -euo pipefail

module purge
module load bcftools
module load vcftools
module load samtools

# ---- inputs/outputs ----
IN_VCF="/lisc/scratch/course/pgbiow/GWAS/Report_DOp25-10208_4.1.vcf"
REF="/lisc/scratch/course/pgbiow/data/genomes/EG5_reference/EG5_reference_genomic.fna"
OUTDIR="/lisc/data/scratch/course/pgbiow/06_diversity_structure"
mkdir -p "$OUTDIR/plink"

CLEAN_PREFIX="$OUTDIR/plink/my_data"   # final cleaned VCF will be ${CLEAN_PREFIX}.vcf.gz

# ---------- Step 0: (optional) choose one replicate per sample ----------
# Keep the lower-missing replicate when IDs are like A_xxx / B_xxx.
# If you always want A_ samples, comment out this whole block and make keep_ids.txt with only A_ IIDs.

vcftools --vcf "$IN_VCF" --missing-indv --out "$OUTDIR/tmp.miss"
awk 'NR>1{
  id=$1; miss=$5; core=id; sub("^[AB]_", "", core);
  if(!(core in best) || miss < best[core]) { best[core]=miss; keep[id]=1 }
}
END{ for (k in keep) print k }' "$OUTDIR/tmp.miss.imiss" > "$OUTDIR/keep_ids.txt"

# Subset to kept samples
vcftools --vcf "$IN_VCF" --keep "$OUTDIR/keep_ids.txt" --recode --stdout | bgzip -c > "$OUTDIR/step0.samples.vcf.gz"
tabix -p vcf "$OUTDIR/step0.samples.vcf.gz"

# ---------- Step 1: split multiallelics, left-align vs reference, drop duplicate records ----------
[ -f "${REF}.fai" ] || samtools faidx "$REF"

bcftools norm -f "$REF" -m -any -d all \
  -Oz -o "$OUTDIR/step1.norm.vcf.gz" "$OUTDIR/step0.samples.vcf.gz"
tabix -p vcf "$OUTDIR/step1.norm.vcf.gz"

# ---------- Step 2: keep only biallelic A/C/G/T SNPs ----------
bcftools view -i 'TYPE="snp" && N_ALT=1 && REF~"^[ACGT]$" && ALT~"^[ACGT]$"' \
  -Oz -o "$OUTDIR/step2.snp_bial.vcf.gz" "$OUTDIR/step1.norm.vcf.gz"
tabix -p vcf "$OUTDIR/step2.snp_bial.vcf.gz"

# ---------- Step 3: apply MAF, missingness, and depth filters with VCFtools ----------
# If your VCF HAS per-genotype FORMAT/DP, use --minDP (per-sample depth):
# If it does NOT, use --min-meanDP (site-average depth) instead.

# Quick check for FORMAT/DP:
HAS_FMT_DP=$(bcftools view -h "$OUTDIR/step2.snp_bial.vcf.gz" | grep -c '##FORMAT=<ID=DP')

if [ "$HAS_FMT_DP" -gt 0 ]; then
  echo "[info] Using per-genotype depth (--minDP)"
  vcftools --gzvcf "$OUTDIR/step2.snp_bial.vcf.gz" \
    --maf 0.05 \
    --max-missing 0.90 \
    --minDP 5 \
    --recode --recode-INFO-all \
    --out "$CLEAN_PREFIX"
else
  echo "[info] FORMAT/DP not found; using site mean depth (--min-meanDP)"
  vcftools --gzvcf "$OUTDIR/step2.snp_bial.vcf.gz" \
    --maf 0.05 \
    --max-missing 0.90 \
    --min-meanDP 5 \
    --recode --recode-INFO-all \
    --out "$CLEAN_PREFIX"
fi

# bgzip + index the final cleaned VCF for downstream convenience
bgzip -c "$CLEAN_PREFIX.recode.vcf" > "$CLEAN_PREFIX.vcf.gz"
tabix -p vcf "$CLEAN_PREFIX.vcf.gz"

echo "Cleaned VCF: $CLEAN_PREFIX.vcf.gz"
```

Submit the job:
```bash
sbatch 00a_clean_vcf.sh
```

### Step 2 — Convert VCF to PLINK format

What this produces:

- `my_data.bed` – binary genotype matrix  
- `my_data.bim` – variant metadata (CHROM, POS, ID, REF/ALT)  
- `my_data.fam` – sample metadata (FID, IID, Sex, Phenotype placeholder)  

**Create the Slurm script with `vi`:**

1. Open the editor:
   ```bash
   vi 00_convert_vcf_to_plink.sh
   ```
2. Press `i` to enter *INSERT* mode.  
3. **Copy & paste** the content below:

```bash
#!/bin/bash
#SBATCH --job-name=vcf2plink_qc
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH -o vcf2plink_qc.out
#SBATCH -e vcf2plink_qc.err

set -euo pipefail

module purge
module load PLINK

# ================= USER SETTINGS =================
# Use the CLEANED VCF produced by 00a_clean_vcf.sh
IN_VCF="/lisc/data/scratch/course/pgbiow/06_diversity_structure/plink/my_data.vcf.gz"
OUT_DIR="/lisc/data/scratch/course/pgbiow/06_diversity_structure/plink"
OUT_BASENAME="my_data"        # final PLINK prefix: ${OUT_DIR}/my_data
MAX_SAMPLE_MISS="0.20"        # drop samples with >20% missing genotypes
THREADS="${SLURM_CPUS_PER_TASK:-1}"
# =================================================

mkdir -p "$OUT_DIR"

echo "[$(date)] Converting cleaned VCF -> PLINK, and filtering high-missing samples ..."
plink --vcf "$IN_VCF" \
      --make-bed \
      --double-id \
      --allow-extra-chr \
      --mind "$MAX_SAMPLE_MISS" \
      --threads "$THREADS" \
      --out "${OUT_DIR}/${OUT_BASENAME}"

echo "[$(date)] Done."
echo "Outputs:"
echo "  - ${OUT_DIR}/${OUT_BASENAME}.bed"
echo "  - ${OUT_DIR}/${OUT_BASENAME}.bim"
echo "  - ${OUT_DIR}/${OUT_BASENAME}.fam"

# Quick sanity prints (optional):
echo -n "Samples: " && wc -l < "${OUT_DIR}/${OUT_BASENAME}.fam"
echo -n "SNPs:    " && wc -l < "${OUT_DIR}/${OUT_BASENAME}.bim"
```
4. Press `Esc`, type `:wq` and press `Enter` to save and exit.  
5. Submit the job:
```bash
sbatch 00_convert_vcf_to_plink.sh
```
---

## Part 1 — Genetic Diversity Estimation (PLINK + VCFtools)

We will compute **heterozygosity**, **inbreeding (F)**, **missingness**, and **allele frequencies** to assess data quality and diversity.  
For consistency, use the **unpruned** dataset (`my_data.*`) for diversity summaries (unless you specifically want summaries on the pruned set).

### Step 1 — PLINK: Heterozygosity and Missingness

Calculating heterozygosity and missingness is important because it allows us to evaluate both the **biological signal** and the **technical quality** of the dataset. Samples with unusually low heterozygosity may indicate **inbreeding** or **loss of diversity**, while samples with extremely high heterozygosity could suggest **contamination** or **sample mix-ups**. Similarly, individuals or loci with high levels of missing data may reflect **low sequencing depth**, **library preparation problems**, or **poor alignment**, and can bias downstream analyses such as PCA, Admixture, or GWAS if not detected and filtered out.

**Create the Slurm script with `vi`:**

1. Open:
   ```bash
   vi 10_plink_diversity.sh
   ```
2. Press `i` to enter *INSERT* mode.  
3. **Copy & paste**:
   ```bash
   #!/bin/bash
   #SBATCH --job-name=plink_div
   #SBATCH --cpus-per-task=1
   #SBATCH --mem=1G
   #SBATCH --time=00:20:00
   #SBATCH -o plink_div.out
   #SBATCH -e plink_div.err

   module load PLINK

   IN_BASE="./plink/my_data"  

   echo "Calculating heterozygosity and inbreeding coefficient (per sample)..."
   plink --bfile "${IN_BASE}" \
         --het \
         --allow-extra-chr \
         --out ./diversity/plink_het

   echo "Calculating per-individual missingness..."
   plink --bfile "${IN_BASE}" \
         --missing \
         --allow-extra-chr \
         --out ./diversity/plink_missing

   echo "Done. Outputs: plink_het.het, plink_missing.imiss"
   ```

4. `Esc`, `:wq`, `Enter`.  
5. Submit:
   ```bash
   sbatch 10_plink_diversity.sh
   ```

you will get two main output files inside the `diversity/` folder:

1) `plink_het.het` — per-individual heterozygosity and inbreeding coefficient
   ```bash
   head diversity/plink_het.het
   ```
Columns:
`O(HOM)` = Observed number of homozygous genotypes
`E(HOM)` = Expected number of homozygous genotypes under Hardy–Weinberg
`N(NM)` = Number of non-missing genotypes
`F` = inbreeding coefficient per individual

**Interpretation of F:**  
**High positive values (~0.1 or higher)** → possible inbreeding or loss of diversity

2) `plink_missing.imiss` — missing data per individual
   ```bash
   head diversity/plink_missing.imiss
   ```
Columns:
`N_MISS` = Number of missing genotypes
`F_MISS` = Fraction of missing data per individual

### Open RStudio on the Cluster
On the LiSC cluster, you do not run RStudio directly from the terminal.
Instead, you:

Open a browser and go to:
https://rstudio.lisc.univie.ac.at

Log in with your cluster username and password.
This gives you an RStudio session running on the cluster.

Once logged in, open a new R script and paste:

```r
# Load heterozygosity results
het <- read.table("diversity/plink_het.het", header=TRUE)

# Load missingness results
miss <- read.table("diversity/plink_missing.imiss", header=TRUE)

# Show plots directly in RStudio
hist(het$F, main="Inbreeding coefficient (F)", xlab="F")
hist(miss$F_MISS, main="Missing genotypes per individual", xlab="Fraction missing")
plot(miss$F_MISS, het$F,
     xlab="Missingness (F_MISS)",
     ylab="Inbreeding coefficient (F)",
     main="QC: Missingness vs. Inbreeding")

# Save all three plots to PDF
pdf("diversity/qc_plots.pdf")

hist(het$F, main="Inbreeding coefficient (F)", xlab="F")
hist(miss$F_MISS, main="Missing genotypes per individual", xlab="Fraction missing")
plot(miss$F_MISS, het$F,
     xlab="Missingness (F_MISS)",
     ylab="Inbreeding coefficient (F)",
     main="QC: Missingness vs. Inbreeding")

dev.off()
```

**Count how many SNPs are actually polymorphic in your data and how informative they are:**  
Create a file called ```panel_qc_stats.sh```and then copy and paste:

```bash
#!/bin/bash
#SBATCH --job-name=panel_qc_stats
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH -o panel_qc_stats.out
#SBATCH -e panel_qc_stats.err

set -euo pipefail

module purge
module load PLINK

# ---- paths ----
BASE="/lisc/data/scratch/course/pgbiow/06_diversity_structure/plink/my_data"   # prefix of .bed/.bim/.fam
OUTD="/lisc/data/scratch/course/pgbiow/06_diversity_structure/diversity"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
mkdir -p "$OUTD"

echo "[$(date)] Step 1: Frequency spectrum"
plink --bfile "$BASE" \
      --allow-extra-chr \
      --freq \
      --threads "$THREADS" \
      --out "$OUTD/panel_freq"

echo "[$(date)] Step 2: Summarize mono/rare/informative"
awk 'NR>1{
  maf=$5;
  mono+=(maf==0);
  rare+=(maf>0 && maf<0.01);
  inf+=(maf>=0.05);
  tot++
}
END{
  printf "Total=%d\nMonomorphic=%d (%.1f%%)\nRare(MAF<0.01)=%d (%.1f%%)\nInformative(MAF>=0.05)=%d (%.1f%%)\n",
         tot, mono, 100*mono/tot, rare, 100*rare/tot, inf, 100*inf/tot
}' "$OUTD/panel_freq.frq" | tee "$OUTD/panel_freq.summary.txt"

echo "[$(date)] Step 3: Recompute F on informative SNPs only (MAF>=0.05)"
plink --bfile "$BASE" \
      --allow-extra-chr \
      --maf 0.05 \
      --make-bed \
      --threads "$THREADS" \
      --out "$OUTD/my_data.maf05"

plink --bfile "$OUTD/my_data.maf05" \
      --allow-extra-chr \
      --het \
      --threads "$THREADS" \
      --out "$OUTD/plink_het_maf05"

echo "[$(date)] Done."
echo "Outputs:"
echo "  - $OUTD/panel_freq.frq"
echo "  - $OUTD/panel_freq.summary.txt"
echo "  - $OUTD/my_data.maf05.{bed,bim,fam}"
echo "  - $OUTD/plink_het_maf05.het"
```

Submit:
```bash
sbatch panel_qc_stats.sh
```

### Post-Run QC (Oil Palm Panel)


**1) Quick sanity: counts**  

```bash
OUTD="/lisc/data/scratch/course/pgbiow/06_diversity_structure/diversity"

# Samples and SNPs after MAF≥0.05 filter
echo "Samples:" $(wc -l < "$OUTD/my_data.maf05.fam")
echo "SNPs:"    $(wc -l < "$OUTD/my_data.maf05.bim")
```

**Interpretation:**
- Make sure the sample count matches expectations.
- SNP count is how many variants survived the **MAF ≥ 0.05** filter.

**2) Frequency spectrum summary (already computed)**  

```bash
# Peek at the first few freq rows
head "$OUTD/panel_freq.frq"

# Monomorphic / rare / informative summary you generated
cat "$OUTD/panel_freq.summary.txt"
```

**Tip:** If “**Informative (MAF≥0.05)**” is low, inflated F in raw data was expected (lots of near-fixed loci).


**3) Inbreeding (F) on informative SNPs**  

```bash
# Preview
head "$OUTD/plink_het_maf05.het"

# Summary stats
awk 'NR>1{n++; F=$6; if(n==1||F<min)min=F; if(n==1||F>max)max=F; sum+=F}
     END{printf "n=%d  meanF=%.4f  minF=%.4f  maxF=%.4f\n", n, sum/n, min, max}'     "$OUTD/plink_het_maf05.het"

# Top 10 highest/lowest F
awk 'NR>1{print $1,$2,$6}' "$OUTD/plink_het_maf05.het" | sort -k3,3nr | head
awk 'NR>1{print $1,$2,$6}' "$OUTD/plink_het_maf05.het" | sort -k3,3n  | head
```

**Guideline (oil palm):** After QC, most **F** should be around **0 to 0.05**. Much higher suggests population structure (Wahlund effect), close relatives, or lingering QC issues—compute F **within groups** and/or remove duplicates/relatives before final interpretation.


### Step 2 — VCFtools: Heterozygosity, Missingness (Optional)

We will use another tool, **VCFtools**, that complements PLINK.  
While PLINK is strong for data preparation and sample-level summaries, VCFtools works directly on the VCF file and provides detailed per-site and per-individual statistics such as MAF, heterozygosity, and missingness. Combining both tools gives us a more complete and reliable picture of genetic diversity and data quality.


**Create the Slurm script with `vi`:**

1. Open:
   ```bash
   vi 11_vcftools_diversity.sh
   ```
2. Press `i`.  
3. **Copy & paste**:
   ```bash
   #!/bin/bash
   #SBATCH --job-name=vcftools_div
   #SBATCH --cpus-per-task=1
   #SBATCH --mem=1G
   #SBATCH --time=00:20:00
   #SBATCH -o vcftools_div.out
   #SBATCH -e vcftools_div.err

   module load vcftools

   IN_VCF=/lisc/scratch/course/pgbiow/GWAS/Report_DOp25-10208_4.1.vcf"  # same input used in Step 0.1

   echo "Per-individual heterozygosity and inbreeding coefficient..."
   vcftools --vcf "${IN_VCF}" \
            --het \
            --out ./diversity/vcftools_het

   echo "Per-individual missingness..."
   vcftools --vcf "${IN_VCF}" \
            --missing-indv \
            --out ./diversity/vcftools_missing

   echo "Done. Outputs: vcftools_het.het, vcftools_missing.imiss"
   ```
4. `Esc`, `:wq`, `Enter`.  
5. Submit:
   ```bash
   sbatch 11_vcftools_diversity.sh
   ```

**Key outputs**

- `vcftools_het.het` → per-sample O(HOM), E(HOM), **F**  
- `vcftools_missing.imiss` → per-sample missingness


### ❓ Diversity QC — Questions to Consider

- Are there individuals with **very high missingness** (e.g., >10–20%)?  
- Do some samples show **unusually low heterozygosity** (possible inbreeding or data issues)?  

> **Relevance:**
> - Estimating heterozygosity, inbreeding, and allele frequencies helps breeders monitor the **health of germplasm collections**.  
> - Detects individuals with **low diversity** that may reduce genetic gain.  
> - Identifies **outliers or contaminated samples** early.  
> - Guides decisions on which individuals to advance or discard in breeding pipelines.

---

## Part 3 — Principal Component Analysis (PCA) with PLINK
Principal Component Analysis (PCA) is a statistical method that reduces the dimensionality of complex datasets while retaining most of the variation. In population genetics, PCA is used to visualize genetic relationships between individuals and identify population structure. Individuals that are genetically similar will cluster together in the PCA plot, while genetically distinct populations will separate.

We will run PCA on the **LD-pruned** dataset produced in **Step 2** (`my_data_pruned.*`).  

**Inputs used here**

- `my_data_pruned.bed/.bim/.fam` (created by `01_ld_pruning.sh`)

**Outputs produced**

- `pca_results.eigenval` — variance explained by each PC  
- `pca_results.eigenvec` — PC scores per individual  

### Step 1 — Run PCA (create script with `vi`)

1. Open:
   ```bash
   vi 20_run_pca.sh
   ```
2. Press `i`.  
3. **Copy & paste**:
   ```bash
   #!/bin/bash
   #SBATCH --job-name=run_pca
   #SBATCH --cpus-per-task=1
   #SBATCH --mem=1G
   #SBATCH --time=00:30:00
   #SBATCH -o run_pca.out
   #SBATCH -e run_pca.err

   module load PLINK

   # LD-pruned dataset from Part 0
   INPUT_BASE_PRUNED="./plink/my_data"

   echo "Running PCA on ${INPUT_BASE_PRUNED}..."
   plink --bfile "${INPUT_BASE_PRUNED}" \
         --pca \
         --allow-extra-chr \
         --out ./pca/pca_results

   echo "Done. Outputs: pca_results.eigenval, pca_results.eigenvec"
   ```
4. `Esc`, `:wq`, `Enter`.  
5. Submit:
   ```bash
   sbatch 20_run_pca.sh
   ```

### PCA Visualization in R

After running PCA with PLINK, we will visualize the results using R.  
Copy the following script into RStudio (or an interactive R session on the cluster) and run it.

*Please change the path to your pca results in the script:


```r
# ================================
# PCA QC for technical replicates (A_/B_)
# ================================

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

# ---- Paths ----
eigvec_path <- "/lisc/scratch/course/pgbiow/06_diversity_structure/pca/pca_results.eigenvec"
eigval_path <- "/lisc/scratch/course/pgbiow/06_diversity_structure/pca/pca_results.eigenval"
out_dir     <- "/lisc/scratch/course/pgbiow/06_diversity_structure/pca"

# ---- Load PCA ----
ev <- read.table(eigvec_path, header = FALSE, stringsAsFactors = FALSE)
colnames(ev) <- c("FID","IID", paste0("PC", seq_len(ncol(ev)-2)))
eigs <- scan(eigval_path, quiet = TRUE)
var_exp <- round(100 * eigs / sum(eigs), 2)

cat("Variance explained (%):\n")
print(var_exp[1:min(10, length(var_exp))])

# ---- Parse replicate info ----
# Replicate label is "A" or "B" if ID starts with "A_" or "B_", else NA
ev <- ev %>%
  mutate(Rep = ifelse(grepl("^[AB]_", IID), substr(IID, 1, 1), NA),
         CoreID = sub("^[AB]_", "", IID))  # remove leading "A_" or "B_"

cat("\nReplicate counts (A/B/NA):\n")
print(table(ev$Rep, useNA = "ifany"))

# ---- Pair up A/B replicates by CoreID for QC ----
A <- ev %>% filter(Rep == "A") %>% select(CoreID, starts_with("PC"))
B <- ev %>% filter(Rep == "B") %>% select(CoreID, starts_with("PC"))

paired <- inner_join(A, B, by = "CoreID", suffix = c("_A","_B"))
cat("\nPaired A/B samples found:", nrow(paired), "\n")

# Compute distances between A and B (PC1..PC10 if available, else up to available PCs)
max_pc <- min(10, sum(grepl("^PC", names(paired))))
pcA <- as.matrix(paired %>% select(all_of(paste0("PC", 1:max_pc, "_A"))))
pcB <- as.matrix(paired %>% select(all_of(paste0("PC", 1:max_pc, "_B"))))
d_pc <- sqrt(rowSums((pcA - pcB)^2))
paired$dist_PC1_2 <- sqrt((paired$PC1_A - paired$PC1_B)^2 + (paired$PC2_A - paired$PC2_B)^2)
paired$dist_PCk   <- d_pc

cat("\nReplicate distance summary (Euclidean):\n")
print(summary(paired$dist_PC1_2))
cat("95th percentile (PC1-2):", quantile(paired$dist_PC1_2, 0.95, na.rm = TRUE), "\n")
cat("95th percentile (PC1-", max_pc, "): ", quantile(paired$dist_PCk, 0.95, na.rm = TRUE), "\n", sep = "")

# Flag outlier pairs (top 5 by PC1-2 distance)
paired_outliers <- paired %>%
  arrange(desc(dist_PC1_2)) %>%
  head(5) %>%
  select(CoreID, dist_PC1_2, dist_PCk)
cat("\nTop 5 A/B pairs by PC1-2 distance:\n")
print(paired_outliers)

# ---- Build a plotting dataframe with points for A and B, and segments connecting pairs ----
plot_df <- ev %>%
  filter(Rep %in% c("A","B")) %>%
  select(IID, CoreID, Rep, PC1, PC2)

segments_df <- paired %>%
  transmute(CoreID,
            xA = PC1_A, yA = PC2_A,
            xB = PC1_B, yB = PC2_B)

# ---- Plot: PC1 vs PC2 with A/B points and connecting segments ----
p <- ggplot() +
  # segments first (so points sit on top)
  geom_segment(data = segments_df,
               aes(x = xA, y = yA, xend = xB, yend = yB),
               linewidth = 0.5, alpha = 0.5) +
  geom_point(data = plot_df,
             aes(x = PC1, y = PC2, color = Rep, shape = Rep),
             size = 2.8, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA with technical replicates (A/B) connected",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

print(p)

ggsave(file.path(out_dir, "pca_replicates_connected_PC1_PC2.pdf"),
       p, width = 8, height = 6)
cat("\nSaved:", file.path(out_dir, "pca_replicates_connected_PC1_PC2.pdf"), "\n")

# ---- Optional: collapse replicates to a single point per CoreID (mean PCs) ----
collapsed <- ev %>%
  group_by(CoreID) %>%
  summarize(across(starts_with("PC"), mean, na.rm = TRUE), .groups = "drop")

p_mean <- ggplot(collapsed, aes(x = PC1, y = PC2)) +
  geom_point(size = 2.6, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA collapsed to per-sample centroids (mean of A/B)",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  )

print(p_mean)

ggsave(file.path(out_dir, "pca_collapsed_centroids_PC1_PC2.pdf"),
       p_mean, width = 8, height = 6)
cat("Saved:", file.path(out_dir, "pca_collapsed_centroids_PC1_PC2.pdf"), "\n")
```

**Interpreting PCA**

- **PC1, PC2** typically separate major population groups; nearby individuals in the plot are genetically similar.  
- Use `pca_results.eigenval` to compute the **variance explained** (e.g., PC1% = eigenval₁ / sum(eigenvals) × 100).  
- The file `pca_results.eigenvec` begins with `FID IID` followed by PC1, PC2, … per individual.

> **Relevance:**  
> - PCA provides a **quick visualization of genetic relationships** among lines.  
> - Helps breeders select **divergent parents** to maximize heterosis.  
> - Detects **substructure** in germplasm panels, preventing bias in association studies.  
> - Highlights **redundant individuals** that could be removed to save costs.
---

## Part 4 - Admixture Analysis

**What is Admixture?**

ADMIXTURE is a software tool for estimating individual ancestries from multilocus SNP genotype datasets. It assumes a model where individuals originate from a mixture of `K` ancestral populations. The output provides the proportion of ancestry each individual derives from each of the `K` ancestral populations. Unlike PCA, which is an unsupervised method, Admixture provides a more direct estimate of ancestral contributions.

We will run ADMIXTURE for different values of `K` (number of ancestral populations) and then determine the optimal `K`.

### Step 0 — Preparing Data for ADMIXTURE

ADMIXTURE requires chromosomes in the `.bim` file to be coded as **plain integers** (`1, 2, 3, …`).  
If your `.bim` file contains labels like `chr1` or `chrLG1`, you must recode them before running ADMIXTURE.  
Otherwise, you will get the error:

**Fix chromosome codes in the `.bim` file**

```bash
cd ./plink

# Backup the original .bim file
cp my_data.bim my_data.bim.backup

# Rewrite the first column: remove "chrLG" prefix and keep only the number
awk '{
  if ($1 ~ /^chrLG[0-9]+$/) {
    gsub("chrLG","",$1)
  }
  print
}' OFS='\t' my_data.bim.backup > my_data.bim
```

### Step 1 - Run ADMIXTURE for different K values

ADMIXTURE takes PLINK BED files as input. We will run it for `K=2` to `K=5` as an example. It is recommended to run ADMIXTURE multiple times for each K value with different random seeds to check for convergence.

First, move back to your folder `06_diversity_structure` using `cd`

**Create a file named `run_admixture.sh` using `vi`:**

```bash
#!/bin/bash
#SBATCH --job-name=admixture_clean_prune_run
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=03:00:00
#SBATCH -o admixture/admixture.out
#SBATCH -e admixture/admixture.err

set -euo pipefail

module purge
module load PLINK
module load admixture

# ----- INPUT: -----
BASE="/path to your user/06_diversity_structure/plink/my_data"
OUTDIR="/path to your user/06_diversity_structure/admixture"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
mkdir -p "$OUTDIR"

echo "[$(date)] Step 0: Remove all-missing and monomorphic loci"
plink --bfile "$BASE" --allow-extra-chr \
      --geno 0.999 --mac 1 \
      --make-bed --out "$OUTDIR/my_data.clean"

echo "[$(date)] Step 1: LD prune on the cleaned set"
plink --bfile "$OUTDIR/my_data.clean" --allow-extra-chr \
      --indep-pairwise 50 5 0.2 \
      --out "$OUTDIR/prune"

echo "[$(date)] Step 2: Build pruned dataset"
plink --bfile "$OUTDIR/my_data.clean" --allow-extra-chr \
      --extract "$OUTDIR/prune.prune.in" \
      --make-bed --out "$OUTDIR/my_data.pruned"

echo "[$(date)] Step 3: Convert contig names to integer chromosome codes (ADMIXTURE requires integers)"
awk '{
  chr=$1; if (!(chr in seen)) { seen[chr]=++c } $1=seen[chr]; print
}' "$OUTDIR/my_data.pruned.bim" > "$OUTDIR/my_data.pruned.intchr.bim"
cp "$OUTDIR/my_data.pruned.bed" "$OUTDIR/my_data.pruned.intchr.bed"
cp "$OUTDIR/my_data.pruned.fam" "$OUTDIR/my_data.pruned.intchr.fam"

# IMPORTANT: ADMIXTURE writes outputs to the *current working directory*
cd "$OUTDIR"

INTBASE="my_data.pruned.intchr"   # basename only; no path

echo "[$(date)] Step 4: Run ADMIXTURE with 5-fold CV (K=2..6)"
for K in 2 3 4 5 6; do
  echo "Running ADMIXTURE K=$K"
  admixture --cv=5 -j"$THREADS" "${INTBASE}.bed" $K | tee "admixture_K${K}.log"

  # Rename outputs to consistent filenames (they are created in $OUTDIR thanks to cd)
  if [[ -f "${INTBASE}.${K}.Q" ]]; then mv "${INTBASE}.${K}.Q" "admixture_K${K}.Q"; fi
  if [[ -f "${INTBASE}.${K}.P" ]]; then mv "${INTBASE}.${K}.P" "admixture_K${K}.P"; fi

  # Show CV error line
  grep -E "CV error" "admixture_K${K}.log" || true
done

echo "[$(date)] Step 5: Summarize CV errors"
grep -H "CV error" admixture_K*.log | tee "CV_errors_summary.txt" || true

echo "[$(date)] Done."
echo "Outputs in $OUTDIR:"
echo "  - my_data.clean.{bed,bim,fam}"
echo "  - my_data.pruned.{bed,bim,fam}"
echo "  - my_data.pruned.intchr.{bed,bim,fam}"
echo "  - admixture_K*.{Q,P,log}"
echo "  - CV_errors_summary.txt"
```

**Submit the job:**

```bash
sbatch run_admixture.sh
```

*   The `--cv` flag performs cross-validation, which helps in determining the optimal `K` value.
*   The `tee` command saves the output to both the screen and a log file.

**❓ Question:** What does the `K` parameter in ADMIXTURE represent? How does the algorithm estimate ancestry proportions?

### Step 2 - Determine Optimal K and Visualize Admixture Results (Conceptual)

To find the optimal `K`, you typically look for the `K` value with the lowest cross-validation error (CV error) from the `.log` files generated by ADMIXTURE.

```bash
# Extract CV errors from log files
grep "CV error" ./admixture/admixture_K*.log
```

### Step 3 — Visualization of ADMIXTURE Results in R

ADMIXTURE results are usually visualized as **stacked bar plots**:  
- Each bar = one individual  
- Each color = ancestry proportion from one of the `K` clusters  
- The height of each color segment = the proportion of ancestry assigned to that cluster  

This gives a clear picture of how individuals are structured genetically and whether populations are admixed.


#### Example R Script

Copy this into your R session (on the cluster with graphics enabled or in RStudio):

```r
# ======================================
# ADMIXTURE barplot (K=6), no metadata
# - orders individuals by major cluster
# - draws cluster blocks + labels
# ======================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ---- Paths (edit if needed) ----
OUTDIR <- "/lisc/data/scratch/course/pgbiow/06_diversity_structure/admixture"
Qfile  <- file.path(OUTDIR, "admixture_K6.Q")                # from your ADMIXTURE run
FAM    <- file.path(OUTDIR, "my_data.pruned.fam")            # same sample order used for ADMIXTURE
OUTPDF <- file.path(OUTDIR, "admixture_K6_sorted_by_cluster.pdf")

# ---- 1) Read Q and FAM ----
Q <- read_table(Qfile, col_names = FALSE, show_col_types = FALSE)
K <- ncol(Q)
colnames(Q) <- paste0("Cluster", seq_len(K))

fam <- read_table(FAM, col_names = FALSE, show_col_types = FALSE)
# FAM columns: FID IID PAT MAT SEX PHENO
colnames(fam) <- c("FID","IID","PAT","MAT","SEX","PHENO")

stopifnot(nrow(Q) == nrow(fam))  # ensure matching sample count

# Attach IDs and compute "major cluster" (argmax)
Q$IID <- fam$IID
Q$FID <- fam$FID
Q <- Q %>%
  relocate(FID, IID)

Q <- Q %>%
  rowwise() %>%
  mutate(
    MajorCluster = paste0("Cluster", which.max(c_across(starts_with("Cluster")))),
    MajorValue   = max(c_across(starts_with("Cluster")))
  ) %>%
  ungroup()

# ---- 2) Order samples by MajorCluster, then by decreasing membership in that cluster ----
ordered_ids <- Q %>%
  arrange(MajorCluster, desc(MajorValue)) %>%
  mutate(Index = row_number()) %>%
  select(FID, IID, MajorCluster, MajorValue, Index)

# Long format for stacked bar plot
Q_long <- Q %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Cluster", values_to = "Ancestry") %>%
  left_join(ordered_ids, by = c("FID","IID"))

# Compute block boundaries and label centers per MajorCluster
block_sizes <- ordered_ids %>%
  count(MajorCluster, name = "n") %>%
  arrange(MajorCluster)  # alphabetical: Cluster1..ClusterK

# In plotting order (which follows arrange(MajorCluster, desc(MajorValue)))
block_sizes <- ordered_ids %>%
  group_by(MajorCluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  # keep the order as it appears in 'ordered_ids'
  mutate(MajorCluster = factor(MajorCluster, levels = unique(ordered_ids$MajorCluster))) %>%
  arrange(MajorCluster)

bounds  <- head(cumsum(block_sizes$n), -1) + 0.5
centers <- cumsum(block_sizes$n) - block_sizes$n/2

# ---- 3) Plot ----
p <- ggplot(Q_long, aes(x = Index, y = Ancestry, fill = Cluster)) +
  geom_col(width = 1) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),     # hide sample labels (too many)
    axis.ticks.x     = element_blank(),
    plot.margin      = margin(t = 10, r = 10, b = 35, l = 10),
    legend.position  = "bottom"
  ) +
  labs(
    title = paste0("ADMIXTURE (K=", K, ") — samples ordered by major cluster"),
    x = "Individuals",
    y = "Ancestry proportion"
  ) +
  # vertical separators between major-cluster blocks
  geom_vline(xintercept = bounds, linetype = "dashed", color = "grey30") +
  # block labels under the bars
  annotate("text", x = centers, y = -0.06,
           label = as.character(block_sizes$MajorCluster),
           vjust = 1, size = 3.2, fontface = "bold") +
  coord_cartesian(ylim = c(-0.1, 1), clip = "off") +
  scale_fill_brewer(palette = "Set2")

print(p)

# ---- 4) Save PDF ----
ggsave(OUTPDF, p, width = 12, height = 5)
cat("Saved plot to:\n", OUTPDF, "\n")

# ---- Optional: write a table with sample order and major cluster ----
write_tsv(ordered_ids, file.path(OUTDIR, "admixture_K6_sample_order.tsv"))
cat("Wrote sample order table:\n", file.path(OUTDIR, "admixture_K6_sample_order.tsv"), "\n")
```

This will generate a **stacked bar plot** of ADMIXTURE proportions for all individuals.


### ❓ Interpretation Questions

1. How would you expect individuals from the same population to look in the ADMIXTURE barplot?  
2. What does it mean if an individual shows **multiple colors** (i.e., ancestry from several clusters)?  
3. If the ADMIXTURE plot for **K=2** shows populations split cleanly, but for **K=3** some groups split further, what might that suggest about population structure?  
4. Why is it important to compare results across multiple values of K rather than relying on a single one?

> **Relevance:**  
> - ADMIXTURE shows the **ancestry composition of individuals**.  
> - Identifies **introgressed material** carrying valuable traits.  
> - Reveals **hidden relatedness** or mislabeling of lines.  
> - Supports **line classification** for better crossing schemes and panel design.
---

## Part 5 - SNP Density Plots (Optional)

**What are SNP Density Plots?**

SNP density plots visualize the distribution of Single Nucleotide Polymorphisms (SNPs) across the genome. They help identify regions with high or low genetic variation, which can be indicative of selective sweeps, recombination hotspots, or genomic rearrangements. These plots are typically generated by counting SNPs within fixed-size genomic windows.

### Step 1 - Prepare Data for SNP Density Plot

We need to extract chromosome and position information from our VCF file and then bin the SNPs into genomic windows. We can use `bcftools query` and some shell scripting for this.

**Create a file named `snp_density_data.sh` using `vi`:**

```bash
#!/bin/bash
#SBATCH --job-name=snp_density_data
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH -o snp_density_data.out
#SBATCH -e snp_density_data.err

module load bcftools

# Input VCF
INPUT_VCF="/lisc/scratch/course/pgbiow/GWAS/Report_DOp25-10208_4.1.vcf"

# Window size (100 kb)
WINDOW_SIZE=100000

# Make sure output folder exists
mkdir -p ./snv_density

echo "Extracting SNP positions and calculating density..."

# Correct bcftools + awk command
bcftools query -f '%CHROM\t%POS\n' "$INPUT_VCF" | \
awk -v WS=$WINDOW_SIZE 'BEGIN {OFS="\t"} {
    win_start = int($2/WS)*WS
    win_end   = win_start + WS - 1
    print $1, win_start, win_end
}' | \
sort -k1,1 -k2,2n | \
uniq -c > ./snv_density/snp_density.txt

echo "SNP density data generated: ./snv_density/snp_density.txt"
```

**Explanation of `awk` command:**

*   `int($2/WS)*WS`: Calculates the start position of the genomic window for the current SNP.
*   `int($2/WS)*WS + WS - 1`: Calculates the end position of the genomic window.
*   `uniq -c`: Counts the number of identical lines (i.e., SNPs falling into the same window).

**Submit the job:**

```bash
sbatch snp_density_data.sh
```

This will create `snp_density.txt`, a file with three columns: `count`, `chromosome`, `window_start`, `window_end`.


### Step 2 — Visualize SNP Density Plot in R

Now that we have `snp_density.txt`, we can plot SNP density across the genome using **R**.

#### Example R Script

Copy this into your R session (on the cluster with graphics enabled or in RStudio):

```r
# ================================
# SNP Density Heatmap
# ================================

library(ggplot2)
library(dplyr)
library(scales)

# 1) Load
snp_df <- read.table("./snv_density/snp_density.txt", header = FALSE)
colnames(snp_df) <- c("Count", "Chromosome", "Window_Start", "Window_End")

# 2) Canonicalize contig names (collapse aliases like "..._chromosome_6_EG5")
snp_df <- snp_df %>%
  mutate(
    Canonical = sub("_chromosome_.*$", "", Chromosome),           # e.g., "NC_025998.1"
    Window_Mid = (Window_Start + Window_End) / 2,
    Window_Width = pmax(1, Window_End - Window_Start)             # avoid zero width
  )

# 3) Identify top 16 contigs by span (max Window_End per Canonical)
top16_tbl <- snp_df %>%
  group_by(Canonical) %>%
  summarise(span = max(Window_End, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(span)) %>%
  slice_head(n = 16)

top16_ids <- top16_tbl$Canonical

snp_top <- snp_df %>%
  filter(Canonical %in% top16_ids)

# 4) Order y-axis by decreasing span (largest at top)
snp_top$Canonical <- factor(
  snp_top$Canonical,
  levels = top16_tbl$Canonical         # already sorted desc
)

# 5) Plot with proper tile width (in Mb)
p <- ggplot(snp_top,
            aes(x = Window_Mid/1e6, y = Canonical, fill = Count)) +
  geom_tile(aes(width = Window_Width/1e6), height = 0.9) +
  scale_fill_gradientn(colours = c("darkgreen","green","yellow","orange","red"),
                       name = "SNPs") +
  scale_x_continuous(labels = label_number(suffix = "Mb")) +
  theme_minimal(base_size = 14) +
  labs(title = "SNP Density per 100 kb window (Top 16 Contigs)",
       x = "Genomic Position",
       y = "Contig") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)

# Save
pdf("snv_density/snp_density_heatmap_top16.pdf", width = 10, height = 6)
print(p)
dev.off()
```

**Output**

- A line plot showing the **number of SNPs per window**, separated by chromosome.  
- Regions with unusually high or low density can highlight **selective sweeps, low-recombination regions, or data artifacts**.

#### ❓ Questions 

1. How does changing the **window size** (e.g., 10 kb vs 500 kb) affect the smoothness of the plot?  
2. What biological processes could explain **regions of low SNP density**?  
3. If one chromosome shows systematically fewer SNPs than others, what could be the cause?  

> **Relevance:**  
> - SNP density highlights **variation-rich and variation-poor regions** across the genome.  
> - Guides **marker design** for breeding (dense SNP coverage where needed).

---

## Conclusion

Congratulations! You have completed a comprehensive practical on population structure analysis and SNP density plotting. You have learned how to prepare VCF data for specialized software, perform PCA and Admixture analysis to infer population relationships and ancestries, and generate SNP density plots to visualize genomic variation. These are powerful techniques for understanding the evolutionary history and genetic diversity of populations.

Remember that these analyses are often iterative, and the interpretation of results requires careful consideration of biological context and potential biases. 

## You have completed **Day 6**!

### Useful Tutorials and Resources

- [PLINK Tutorial (Population Genetics)](https://zzz.bwh.harvard.edu/plink/tutorial.shtml)  
- [PCA and Population Structure](https://speciationgenomics.github.io/pca/)  
- [PCA in R](https://rpubs.com/nanette/PCA)  
- [VCFtools User Guide](https://vcftools.sourceforge.net/man_latest.html)
- [Population Genetics - ADMIXTURE](https://speciationgenomics.github.io/ADMIXTURE/)
