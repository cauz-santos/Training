## Week 2: Applied Bioinformatics for Genomics and Breeding

## Day 7: SNP Applications in Breeding - GWAS and Gene Interpretation

Welcome to Day 7!

Today we connect **SNP genotypes** to a **binary** and **quantitative trait** (plant disease infection response) using **Genome-Wide Association Studies (GWAS)**, and then interpret significant signals biologically by mapping top SNPs to **genes** and **functions**.

This is a **hands-on** session; you will:
- Prepare genotype/phenotype/covariates for GWAS
- Run GWAS in **PLINK** (with population-structure correction via **PC covariates**)
- Visualize results (**Manhattan** and **QQ** plots in R)
- Identify **top SNPs**, map them to **genes** (GFF3), and extract **putative functions**
- Discuss how hits feed into **Marker-Assisted Selection (MAS)** and **Genomic Selection (GS)**

> **Relevance:**  
> GWAS uncovers **genomic regions** that control traits like plant disease response. These become **markers** for MAS, and inform **genome-wide prediction** in GS—shortening cycles and improving accuracy in breeding pipelines.

---

### Part 1 - Prepare data

### Step 0 — Prepare Genotypes, Phenotypes, and Covariates

First create a folder for the file outputs of day 7 in your home directory:
   ```bash
   mkdir 07_gwas_selection
   ```
Now, enter the folder with the command `cd`

Then create some subfolders, for each specific analysis we will perform:
   ```bash
   mkdir plink gwas mas_markers gs-lite
   ```

### Step 1 — Build a PLINK-friendly phenotype file for **infection**

We’ll construct two phenotype files:  
- `pheno_infected.txt` → for binary trait
- `pheno_audpc.txt` → for quantitative trait

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

PHENO_CSV="/lisc/scratch/course/pgbiow/GWAS/phenotypes_variant_verdant.csv"  # columns: No.,ID,Infected_Status,AUDPC,Internal_Symptoms
FAM="/lisc/scratch/course/pgbiow/GWAS/plink/data_pruned.fam"

# Ensure stable sort
export LC_ALL=C

############################
# Infected_Status (binary) #
# PLINK expects: 1 = control, 2 = case, -9 = missing
############################
echo -e "FID\tIID\tINFECTED_STATUS" > pheno_infected_12.txt

# Map CSV -> (IID, code) with 1/2/-9 encoding
awk -F',' 'NR>1 {
  gsub(/\r$/,"")                        # strip CR if present
  id=$2; v=$3
  if (v=="" || v=="NA")      code=-9
  else if (v==0)             code=1     # control
  else if (v==1)             code=2     # case
  else                       code=-9
  print id "\t" code
}' "${PHENO_CSV}" | sort -k1,1 > infected_map.tsv

# Join to FAM order; FID=first col, IID=second col in .fam
awk 'NR==FNR { ph[$1]=$2; next }
     { fid=$1; iid=$2; val=(iid in ph ? ph[iid] : -9);
       print fid "\t" iid "\t" val }' \
    infected_map.tsv "${FAM}" >> pheno_infected_12.txt

rm -f infected_map.tsv

################
# AUDPC (quant) #
# Use -9 for missing/non-numeric
################
echo -e "FID\tIID\tAUDPC" > pheno_audpc.txt

# Map CSV -> (IID, value), sanitize to numeric; else -9
awk -F',' 'NR>1 {
  gsub(/\r$/,"")
  id=$2; v=$4
  if (v=="" || v=="NA") { val="-9" }
  else if (v ~ /^-?[0-9]+(\.[0-9]+)?$/) { val=v }  # numeric
  else { val="-9" }
  print id "\t" val
}' "${PHENO_CSV}" | sort -k1,1 > audpc_map.tsv

awk 'NR==FNR { ph[$1]=$2; next }
     { fid=$1; iid=$2; val=(iid in ph ? ph[iid] : -9);
       print fid "\t" iid "\t" val }' \
    audpc_map.tsv "${FAM}" >> pheno_audpc.txt

rm -f audpc_map.tsv

################
# Summaries
################
echo "Done:"
echo "  - pheno_infected_12.txt (FID IID INFECTED_STATUS; 1=control, 2=case, -9=missing)"
echo "  - pheno_audpc.txt       (FID IID AUDPC; -9=missing)"

# Quick counts (optional)
echo "Counts (infected_12):"
tail -n +2 pheno_infected_12.txt | awk '{c[$3]++} END{for(k in c) print k, c[k]}' | sort -k1,1
echo "Missing in AUDPC:"
tail -n +2 pheno_audpc.txt | awk '$3==-9{m++} END{print (m?m:0)}'
```

```bash
sbatch 01_make_pheno.sh
```

### Step 2 — Prepare **PCA covariates** (PC1–PC10) 
Add PC1–PC5 as covariates so GWAS controls for broad genetic background (ancestry/relatedness) differences among samples. Without them, some SNPs look “significant” just because certain groups both share those alleles and tend to have different sucrose values. PCs capture that group effect; adjusting for them removes this bias and leaves signals that are more likely to be truly linked to sucrose.

Please from the directory of day 7 enter the directory `plink`:

```bash
cd plink
```

Create **`covar_pcs.txt`** with header `FID IID PC1 PC2 PC3 PC4 PC5 ... PC10`:

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

set -euo pipefail

EIGENVEC="/lisc/scratch/course/pgbiow/GWAS/plink/cov_pca.eigenvec"   # FID IID PC1 PC2 ...
NUM_PCS=10                                                           # <-- change to 20 if you like
OUTFILE="covar_pcs${NUM_PCS}.txt"

# Build header + rows: FID IID PC1..PCk
awk -v k="${NUM_PCS}" 'BEGIN{
  OFS="\t";
  printf "FID\tIID";
  for(i=1;i<=k;i++) printf "\tPC%d", i;
  printf "\n";
}
{
  OFS="\t";
  out = $1 OFS $2;
  for(i=1;i<=k;i++) out = out OFS $(2+i);
  print out;
}' "${EIGENVEC}" > "${OUTFILE}"

echo "Done: ${OUTFILE} (PC1..PC${NUM_PCS})"
```

```bash
sbatch 02_make_covariates.sh
```

> **Why PCs as covariates?**  
> PCs capture **population structure**; including them reduces **false positives** (spurious associations due to stratification).

---

## Copying PLINK Files to Your User Folder

To copy all necessary PLINK files from the course directory to your own PLINK working directory, use the command below:

```bash
cp /lisc/scratch/course/pgbiow/GWAS/plink/* ./plink/
```

### Part 2 — GWAS for Disease Traits with PLINK (logistic and linear models)

We’ll run:
A **logistic regression** for the **binary trait** Infected_Status
A **linear regression** for the **quantitative trait** AUDPC

Both models include PC1–PC5 as covariates to control for population structure (ancestry/relatedness), which reduces false positives and helps detect true genotype–phenotype associations.


### Which GWAS Regression to Use?

Choose the regression model based on the type of your trait:

#### 1. Linear Regression — `--linear`

Use this when your trait is **quantitative** (i.e., numerical and continuous), such as measurements of disease severity, height, yield, or enzyme activity.  
It estimates how the trait value changes depending on the number of copies of a genetic variant.

#### 2. Logistic Regression — `--logistic`

Use this when your trait is **binary**, with two categories such as infected vs. uninfected or resistant vs. susceptible.  
It estimates how the probability of being in one category (e.g., infected) changes depending on the genetic variant.

#### Summary Table

| Trait            | Type          | PLINK Option | Output File           | Model    |
|------------------|---------------|--------------|------------------------|----------|
| `Infected_Status`| Binary (0/1)  | `--logistic` | `*.assoc.logistic`     | Logistic |
| `AUDPC`          | Quantitative  | `--linear`   | `*.assoc.linear`       | Linear   |


**Step 1 - Infected_Status → Logistic Regression**  

Create the slurm file:
```bash
vi 20_run_gwas_infected.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=gwas_infected
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH -o gwas_infected.out
#SBATCH -e gwas_infected.err

set -euo pipefail

module load PLINK
mkdir -p ./gwas

plink \
  --bfile ./plink/data_pruned \
  --pheno ./pheno_infected_12.txt \
  --covar ./plink/covar_pcs10.txt \
  --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --logistic hide-covar --ci 0.95 \
  --allow-no-sex \
  --allow-extra-chr \
  --threads 8 \
  --out ./gwas/gwas_infected_pc10

echo "Done: gwas_infected_pc10.assoc.logistic"
```

```bash
sbatch 20_run_gwas_infected.sh
```

**Step 2 - AUDPC → Linear Regression**  
Create the slurm file:
```bash
vi 21_run_gwas_audpc.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=gwas_audpc
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH -o gwas_audpc.out
#SBATCH -e gwas_audpc.err

set -euo pipefail

module load PLINK
mkdir -p ./gwas

plink \
  --bfile ./plink/data_pruned \
  --pheno ./pheno_audpc.txt \
  --covar ./plink/covar_pcs10.txt \
  --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --linear hide-covar \
  --allow-no-sex \
  --allow-extra-chr \
  --threads 8 \
  --out ./gwas/gwas_audpc_pc10

echo "Done: gwas_audpc_pc10.assoc.linear"
```

```bash
sbatch 21_run_gwas_audpc.sh
```


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
# GWAS Visualization: (Manhattan plot)
# ================================

setwd("/path to your home directory/07_gwas_selection/")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

make_gwas_plots_top16 <- function(infile, trait_label, outdir = "./gwas") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  dt <- fread(infile)
  if ("TEST" %in% names(dt)) dt <- dt[TEST == "ADD"]

  stopifnot(all(c("CHR","BP","P") %in% names(dt)))
  dt[, P := as.numeric(P)]
  dt <- dt[!is.na(P) & !is.na(BP)]

  # Numeric chromosome index (e.g., 'chr01' -> 1)
  dt[, CHR_RAW := as.character(CHR)]
  dt[, CHR_NUM := suppressWarnings(as.numeric(gsub("[^0-9]+", "", CHR_RAW)))]
  if (all(is.na(dt$CHR_NUM))) dt[, CHR_NUM := suppressWarnings(as.numeric(as.character(CHR)))]
  dt <- dt[!is.na(CHR_NUM)]

  # Keep exactly the 16 largest chromosomes by span
  chr_stats <- dt[, .(min_bp = min(BP), max_bp = max(BP)), by = .(CHR, CHR_NUM)]
  chr_stats[, span := pmax(0, max_bp - min_bp)]
  keep_chr <- chr_stats[order(-span)][1:min(16, .N), CHR]
  kept <- dt[CHR %in% keep_chr]
  dropped <- dt[!CHR %in% keep_chr]
  message(sprintf("Keeping %d chromosomes (top by span), dropping %d.",
                  uniqueN(kept$CHR), uniqueN(dropped$CHR)))

  # Lambda and thresholds
  setorder(kept, CHR_NUM, BP)
  kept[, CHISQ := qchisq(1 - P, df = 1)]
  lambda <- median(kept$CHISQ, na.rm = TRUE) / 0.456
  M <- nrow(kept)
  bonf <- 0.05 / max(M, 1)

  # Cumulative genomic position for Manhattan
  chr_sizes <- kept[, .(chr_len = max(BP, na.rm = TRUE)), by = .(CHR, CHR_NUM)]
  setorder(chr_sizes, CHR_NUM)
  chr_sizes[, chr_start := cumsum(shift(chr_len, fill = 0))]
  kept <- merge(kept, chr_sizes[, .(CHR, chr_start)], by = "CHR", all.x = TRUE)
  kept[, BPcum := BP + chr_start]

  axis_df <- kept[, .(center = (min(BPcum) + max(BPcum)) / 2), by = .(CHR, CHR_NUM)][order(CHR_NUM)]
  kept[, logP := -log10(P)]
  kept[, color := as.factor(CHR_NUM %% 2)]

  # Manhattan
  p_manhattan <- ggplot(kept, aes(BPcum, logP, color = color)) +
    geom_point(alpha = 0.75, size = 0.8) +
    geom_hline(yintercept = -log10(bonf), linetype = "dashed", color = "red", linewidth = 0.6) +
    scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR_NUM) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Chromosome", y = "-log10(p)",
         title = paste("GWAS Manhattan —", trait_label, "(top 16 chromosomes)")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  # QQ
  observed <- sort(kept$P)
  expected <- -log10(ppoints(length(observed)))
  observed_logp <- -log10(observed)
  p_qq <- ggplot(data.frame(Expected = expected, Observed = observed_logp),
                 aes(Expected, Observed)) +
    geom_point(size = 0.9, alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = sprintf("GWAS QQ — %s (lambda = %.3f)", trait_label, lambda),
         x = "Expected -log10(p)", y = "Observed -log10(p)") +
    theme_minimal(base_size = 13)

  # Show in RStudio
  print(p_manhattan)
  print(p_qq)
  message(sprintf("Lambda: %.3f   Bonferroni: %.3g   SNPs plotted: %d", lambda, bonf, M))

  # ---- Save plots to PDF ----
  ggsave(file.path(outdir, paste0("GWAS_", trait_label, "_Manhattan.pdf")),
         plot = p_manhattan, width = 10, height = 5)
  ggsave(file.path(outdir, paste0("GWAS_", trait_label, "_QQ.pdf")),
         plot = p_qq, width = 6, height = 6)

  # ---- Export top hits ----
  kept_sorted <- copy(kept)[order(P)]
  fwrite(kept_sorted[1:min(20, .N)],
         file.path(outdir, paste0("top20_hits_", trait_label, ".tsv")), sep = "\t")
  fwrite(kept_sorted[P <= bonf],
         file.path(outdir, paste0("bonferroni_hits_", trait_label, ".tsv")), sep = "\t")

  cat(sprintf("PDFs saved: GWAS_%s_Manhattan.pdf, GWAS_%s_QQ.pdf\n", trait_label, trait_label))
  cat(sprintf("Tables saved: top20_hits_%s.tsv, bonferroni_hits_%s.tsv\n", trait_label, trait_label))

  invisible(list(kept_chr = unique(kept$CHR), dropped_chr = unique(dropped$CHR),
                 manhattan = p_manhattan, qq = p_qq,
                 lambda = lambda, bonf = bonf))
}

# Run interactively for each trait:
make_gwas_plots_top16("./gwas/gwas_audpc_pc10.assoc.linear", "AUDPC")
make_gwas_plots_top16("./gwas/gwas_infected_pc10.assoc.logistic", "Infected_Status")
```

#### QQ plot interpretation

- **AUDPC (λ = 0.932):** Slight deflation (conservative p-values). The bulk of points follows the null line; the right-hand tail shows putative true associations. Safe calibration for teaching analyses.
- **Infected_Status (λ = 1.052):** Well-calibrated; points track the null with a clear tail of significant SNPs. Residual stratification appears to be controlled after adding PC1–PC10.

**Notes:** We restricted to the 16 largest chromosomes to avoid scaffold noise. Bonferroni thresholds were ~4.6–4.7×10⁻⁷ given ~106–108k SNPs tested.


### Inspecting GWAS Bonferroni Hits (AUDPC trait)

After running the analyses, TSV files with SNPs passing Bonferroni (α = 0.05) are saved in `./gwas`:

- `bonferroni_hits_AUDPC.tsv` (from linear GWAS)
- `bonferroni_hits_Infected_Status.tsv` (from logistic GWAS)

### Quick look

```bash
cd /path to your folder/07_gwas_selection/gwas

# AUDPC (quantitative)
head bonferroni_hits_AUDPC.tsv
cut -f1,2,3,9 bonferroni_hits_AUDPC.tsv | head   # CHR BP SNP P
```

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
# ================================
# Regional Association Plot (±250 kb) around the Top AUDPC Hit
# ================================
library(data.table)
library(ggplot2)
library(ggrepel)

regional_plot_tophit <- function(infile, trait_label = "AUDPC", outdir = "./gwas", flank_kb = 250) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  dt <- fread(infile)

  # Keep additive model only (PLINK outputs)
  if ("TEST" %in% names(dt)) dt <- dt[TEST == "ADD"]
  stopifnot(all(c("CHR","BP","P") %in% names(dt)))

  dt[, P := as.numeric(P)]
  dt <- dt[!is.na(P) & !is.na(BP) & !is.na(CHR)]
  setorder(dt, P)
  lead <- dt[1]
  chr  <- as.character(lead$CHR)
  pos  <- as.integer(lead$BP)

  message(sprintf("Top SNP: CHR=%s  BP=%s  P=%.3e", chr, format(pos, big.mark=","), lead$P))

  flank <- as.integer(flank_kb) * 1000L
  win <- dt[CHR == chr & BP >= (pos - flank) & BP <= (pos + flank)]
  if (nrow(win) < 50) {  # auto-expand if too few points
    flank <- 500000L
    win <- dt[CHR == chr & BP >= (pos - flank) & BP <= (pos + flank)]
    message(sprintf("Few points in window; expanded to ±%d kb", flank/1000L))
  }

  top_win <- win[which.min(P)]

  p_regional <- ggplot(win, aes(x = BP, y = -log10(P))) +
    geom_point(alpha = 0.8, size = 1) +
    geom_vline(xintercept = pos, linetype = "dashed", color = "red") +
    geom_text_repel(
      data = top_win,
      aes(label = ifelse(SNP == ".", paste0(chr, ":", BP), SNP)),
      size = 3, box.padding = 0.3, segment.color = "grey50", max.overlaps = Inf
    ) +
    labs(
      title = sprintf("Regional Association: %s - %s ±%d kb", trait_label, chr, flank/1000L),
      subtitle = sprintf("Lead at %s:%s   p=%.2e", chr, format(pos, big.mark=","), top_win$P),
      x = "Genomic Position (bp)", y = "-log10(p)"
    ) +
    theme_minimal(base_size = 14)

  print(p_regional)

  base <- file.path(outdir, sprintf("GWAS_%s_Regional_%s_%s",
                                    trait_label, gsub("[^A-Za-z0-9_.-]","_", chr), pos))
  ggsave(paste0(base, ".png"), plot = p_regional, width = 10, height = 4, dpi = 150)
  ggsave(paste0(base, ".pdf"), plot = p_regional, width = 10, height = 4)
  message(sprintf("Saved: %s.png and %s.pdf", base, base))
  invisible(p_regional)
}

# Run ONLY for AUDPC (use your pc10 result file)
regional_plot_tophit("./gwas/gwas_audpc_pc10.assoc.linear", "AUDPC")
```
---

### Part 4 — From SNP to Gene (annotation & function)

Goal: For our **top SNP(s)**, find the **overlapping/nearest gene(s)** and a **putative function**.

> **You need:**  
> - The **same reference genome** used for alignment & variant calling (FASTA)  
> - Its **annotation** in **GFF3** (gene coordinates + attributes, e.g., `ID`, `Name`, `product`)  
> - **BEDTools** module (for genomic intersections)

#### Files for Oil Palm Genome - EG5

These files are available on the cluster:  
`/lisc/scratch/course/pgbiow/data/genomes/EG5_reference`  
EG5_reference_genomic.fna # reference genome (FASTA)  
EG5_genomic.gff  # annotation in GFF3 format  
 
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
#SBATCH --job-name=annot_EG5
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH -o annot_EG5.out
#SBATCH -e annot_EG5.err

set -euo pipefail
export LC_ALL=C

# --- Paths (AUDPC only) ---
GENOME_GFF="/lisc/scratch/course/pgbiow/data/genomes/EG5_reference/EG5_genomic.gff"
OUTDIR="annotation_EG5"
mkdir -p "$OUTDIR"

# Prefer Bonferroni hits; if empty/missing, use top20
BONF="./gwas/bonferroni_hits_AUDPC.tsv"
TOP20="./gwas/top20_hits_AUDPC.tsv"

if [[ -s "$BONF" ]]; then
  IN_GWAS="$BONF"
  echo "Using Bonferroni-significant AUDPC hits: $IN_GWAS"
elif [[ -s "$TOP20" ]]; then
  IN_GWAS="$TOP20"
  echo "Bonferroni file empty/missing; using top-20 AUDPC hits: $IN_GWAS"
else
  echo "Error: Neither $BONF nor $TOP20 exists or is non-empty." >&2
  exit 1
fi

# --- 1) SNPs → BED (keep NC_* accessions exactly; 0-based start) ---
# Expect columns include CHR and BP (from your TSVs written by the plotting helper)
# Header line starts with 'CHR'
awk 'BEGIN{FS=OFS="\t"} NR>1 {
  chr=$1; bp=$3; snp=$2;
  if (bp ~ /^[0-9]+$/) print chr, bp-1, bp, snp;
}' "$IN_GWAS" > "$OUTDIR/top_snps.bed"

# --- 2) Genes (GFF3) → BED using EG5 seqids (e.g., NC_025995.1) ---
# Grab only 'gene' features; build a compact name field id|name|product (best-effort).
awk -F'\t' 'BEGIN{OFS="\t"}
  $3=="gene" {
    chr=$1; start=$4-1; end=$5; attrs=$9;
    id="NA"; name="NA"; prod="NA"; note="NA";
    n=split(attrs,a,";");
    for(i=1;i<=n;i++){
      split(a[i],kv,"=");
      if(kv[1]=="ID" && kv[2]!="") id=kv[2];
      if(kv[1]=="Name" && kv[2]!="") name=kv[2];
      if(kv[1]=="gene" && name=="NA" && kv[2]!="") name=kv[2];
      if(kv[1]=="product" && kv[2]!="") prod=kv[2];
      if(kv[1]=="Note" && kv[2]!="") note=kv[2];
      if(kv[1]=="description" && kv[2]!="") note=kv[2];
    }
    if(prod=="NA" && note!="NA") prod=note;
    print chr, start, end, id "|" name "|" prod;
  }' "$GENOME_GFF" > "$OUTDIR/genes.bed"

# --- 3) Sort & intersect/closest with BEDTools ---
sort -k1,1 -k2,2n "$OUTDIR/top_snps.bed"  > "$OUTDIR/top_snps.sorted.bed"
sort -k1,1 -k2,2n "$OUTDIR/genes.bed"     > "$OUTDIR/genes.sorted.bed"

module load BEDTools

# Overlaps (SNP falls inside a gene)
bedtools intersect -a "$OUTDIR/top_snps.sorted.bed" -b "$OUTDIR/genes.sorted.bed" -wa -wb \
  > "$OUTDIR/snps_in_genes.tsv"

# Nearest gene to each SNP (report first tie)
bedtools closest -a "$OUTDIR/top_snps.sorted.bed" -b "$OUTDIR/genes.sorted.bed" -d -t first \
  > "$OUTDIR/snps_nearest_genes.tsv"

echo "Wrote:"
echo "  $OUTDIR/top_snps.bed"
echo "  $OUTDIR/genes.bed"
echo "  $OUTDIR/snps_in_genes.tsv"
echo "  $OUTDIR/snps_nearest_genes.tsv"
```

```bash
sbatch 30_prepare_bed_for_annotation.sh
```

**Interpreting outputs:**
- `snps_in_genes_with_product.tsv` → if a SNP lies **within** a gene (exon/intron span)  
- `snps_nearest_genes_with_product.tsv` → the **closest** gene and the **distance** (0 if overlapping)

### SNP Annotation Results
We have annotated the top GWAS SNPs with nearby genes using the date palm genome annotation.

Two result files are available inside the `annotation/` folder:

**1) SNPs located **inside genes****  

```bash
cat annotation_EG5/snps_in_genes_with_product.tsv
```

Each line shows:  
- SNP coordinates and ID
- The gene it overlaps
- Gene name and product (function, if available)

This file tells us which SNPs directly fall within annotated genes.

**2) The nearest gene to each SNP**  
```bash
cat annotation_EG5/snps_nearest_genes_with_product.tsv
```

Go to the Uniprot database and look for the potential function of these genes: 

https://www.uniprot.org/

Search for:
```bash
LOC105042204 Elaeis guineensis
LOC105042195 Elaeis guineensis
```

Each line shows:  
- SNP coordinates and ID
- The closest gene
- Gene name and product (function, if available)
- Distance in base pairs

This file helps identify candidate genes near but not directly overlapping the SNP.

**Interpretation:**  
This top SNP sits near HSFA2c, a transcription factor that boosts HSP70/HSP90 chaperones. Those chaperones fold key immune receptors, so small changes in HSFA2 activity can indirectly tune disease resistance, especially under heat, when both the pathogen and host heat-responses are active.

> **If function is missing:**  
> - Your GFF annotation may not include `product`—try `Name`/`gene` attributes or consult the genome’s annotation README.  
> - For deeper function (domains/GO), run tools like **InterProScan** offline on the protein FASTA of candidate genes (advanced, not covered here).

> **Relevance:**  
> - Translating SNPs to **genes and functions** turns statistical signals into **biological hypotheses**, which alleles/genes to track, validate, and deploy in breeding.
> - Accelerates **marker development** and **IP differentiation** around candidate genes/alleles.
> - Informs **assay design** (haplotype tagging, primers) and prioritizes **functional validation** (expression, knockouts).
> - Creates a pipeline from **statistical hit → deployable marker → breeding action**.


### Optional Exercise: Map GWAS Peaks to Nearby Genes (±50 kb, Oil Palm EG5)

#### Goal
Take the lead **AUDPC** SNPs and list **all genes within ±50 kb** of each peak on the **oil palm EG5** reference. This bridges statistical association (SNPs) to biological candidates (genes/functions).

#### Why ±50 kb?
It’s a compact teaching window that usually captures the nearest causal candidates without dragging in too many distant genes. (You can tighten to ±25 kb or expand later; LD-based windows are an advanced follow-up.)

```bash
vi genes_50kb.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=genes_50kb
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:02:00
#SBATCH -o genes_50kb.out
#SBATCH -e genes_50kb.err

set -euo pipefail
export LC_ALL=C

# Paths
GENOME_GFF="/lisc/scratch/course/pgbiow/data/genomes/EG5_reference/EG5_genomic.gff"
OUTDIR="annotation_EG5"
mkdir -p "${OUTDIR}"

# --- 1) Define the SNP peaks (CHR  POS); 50 kb window each
cat > "${OUTDIR}/snps_peaks.tsv" <<EOF
NC_025995.1	49444165
NC_025995.1	49444188
EOF

# --- 2) Make gene BED from EG5 GFF (if missing)
if [[ ! -s "${OUTDIR}/genes.sorted.bed" ]]; then
  awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="gene"{
      chr=$1; start=$4-1; end=$5; attrs=$9;
      id="NA"; name="NA"; prod="NA"; note="NA";
      n=split(attrs,a,";");
      for(i=1;i<=n;i++){
        split(a[i],kv,"=");
        if(kv[1]=="ID" && kv[2]!="") id=kv[2];
        if(kv[1]=="Name" && kv[2]!="") name=kv[2];
        if(kv[1]=="gene" && name=="NA" && kv[2]!="") name=kv[2];
        if(kv[1]=="product" && kv[2]!="") prod=kv[2];
        if(kv[1]=="Note" && kv[2]!="") note=kv[2];
        if(kv[1]=="description" && kv[2]!="") note=kv[2];
      }
      if(prod=="NA" && note!="NA") prod=note;
      print chr, start, end, id "|" name "|" prod;
    }' "${GENOME_GFF}" | sort -k1,1 -k2,2n > "${OUTDIR}/genes.sorted.bed"
fi

# --- 3) Build 50 kb windows around SNPs; merge overlapping windows
awk 'BEGIN{OFS="\t"}
     {chr=$1; pos=$2; start=pos-50000; if(start<0) start=0; end=pos+50000; print chr,start,end,chr":"pos}' \
    "${OUTDIR}/snps_peaks.tsv" \
  | sort -k1,1 -k2,2n > "${OUTDIR}/windows_50kb.bed"

module load BEDTools
bedtools merge -i "${OUTDIR}/windows_50kb.bed" -c 4 -o collapse > "${OUTDIR}/windows_50kb.merged.bed"

# --- 4) Intersect windows with genes → pretty table
bedtools intersect -a "${OUTDIR}/windows_50kb.merged.bed" -b "${OUTDIR}/genes.sorted.bed" -wb \
  > "${OUTDIR}/genes_in_50kb_raw.tsv"

# Format: region_chr  region_start  region_end  snps_in_window  gene_chr  gene_start  gene_end  gene_id  gene_name  product
awk 'BEGIN{OFS="\t"; print "window_chr","window_start","window_end","snps","gene_chr","gene_start","gene_end","gene_id","gene_name","product"}
     {
       split($4, s, ","); snps=$4;                 # collapsed SNP list for the window
       split($8, f, "|"); gid=f[1]; gname=f[2]; prod=f[3];
       print $1, $2, $3, snps, $5, $6, $7, gid, gname, prod;
     }' "${OUTDIR}/genes_in_50kb_raw.tsv" \
  | sort -k5,5 -k6,6n -k8,8 \
  | uniq \
  > "${OUTDIR}/genes_in_50kb_around_peaks.tsv"

echo "Wrote ${OUTDIR}/genes_in_50kb_around_peaks.tsv"
```

Submit the job:
```bash
sbatch genes_50kb.sh
```

Check the output results:
```bash
column -t -s $'\t' annotation_EG5/genes_in_50kb_around_peaks.tsv | less -S
```

---

### Part 5 — From GWAS to Selection Decisions

### Quick GS: Compute a Polygenic Score (PGS) for Oil Palm (EG5) — AUDPC (Infection, Quantitative Trait)  
After GWAS, we know which SNPs are significantly associated with our trait. In **Marker-Assisted Selection (MAS)** we could take a few diagnostic SNPs. But in **Genomic Selection (GS)** we want to use many SNPs together to rank lines.  

Here we build a **Polygenic Score (PGS)** for **AUDPC**: for each plant, we multiply its genotype (0/1/2 copies of the effect allele) by the SNP’s estimated effect size (**BETA**) from the AUDPC GWAS, and sum across selected SNPs. This yields one number per line—the **genomic score**. If the PGS correlates with observed **AUDPC**, we can prioritize candidates accordingly. **Because lower AUDPC means less disease, lines with *lower* PGS are preferred** (when BETA was fit to raw AUDPC and positive BETAs increase AUDPC). If you transformed AUDPC or defined the effect differently, interpret the sign accordingly.


####  Prerequisites

- Genotypes (PLINK bed/bim/fam) after QC and pruning:
  - `./plink/data_pruned.{bed,bim,fam}` (oil palm; chromosome names like `NC_025995.1`)
- Phenotype (AUDPC):
  - `./pheno_audpc.txt` with **3 columns**: `FID IID AUDPC`
- GWAS results (linear model with PC1–PC10):
  - `./gwas/gwas_audpc_pc10.assoc.linear`
  - Helper tables produced earlier: 
    - `./gwas/bonferroni_hits_AUDPC.tsv` (preferred) **or** 
    - `./gwas/top20_hits_AUDPC.tsv` (fallback)

> We assume you’re working in `/path to your home directory/07_gwas_selection`. Adjust paths if needed.


### Step 1 — Run GS-lite script (build weights, compute PGS, join with AUDPC)

Create the file and submit **`50_gs_lite_audpc.sh`**:

```bash
vi 50_gs_lite_audpc.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=gs_lite_audpc
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH -o gs_lite_audpc.out
#SBATCH -e gs_lite_audpc.err

set -euo pipefail

OUTDIR="gs-lite"
GWAS_HITS="./gwas/bonferroni_hits_AUDPC.tsv"   # from your GWAS results
BFILE="plink/data_pruned"                      # PLINK prefix (bed/bim/fam)
PHENO="pheno_audpc.txt"                        # FID IID AUDPC

mkdir -p "$OUTDIR"

echo "[GS-lite] Using Bonferroni hits for weights: $GWAS_HITS"

# 1) Build PLINK --score weights: SNP  A1(effect)  BETA
awk 'BEGIN{FS=OFS="\t"} NR>1 {print $2,$4,$7}' "$GWAS_HITS" > "$OUTDIR/gs_weights.tsv"
echo "[GS-lite] Wrote: $OUTDIR/gs_weights.tsv (SNP  A1  BETA)"

# 2) Compute PGS with PLINK
module load PLINK
plink --bfile "$BFILE" \
      --allow-extra-chr \
      --score "$OUTDIR/gs_weights.tsv" 1 2 3 header sum \
      --out "$OUTDIR/gs_pgs"

# 3) Rank and correlate with phenotype in R
module load R
R --vanilla <<'EOF'
prof <- read.table("gs-lite/gs_pgs.profile", header=TRUE, stringsAsFactors=FALSE)
score_col <- grep("^SCORE", names(prof), value=TRUE)[1]
pgs <- prof[, c("FID","IID", score_col)]
names(pgs) <- c("FID","IID","PGS")
pgs$PGS <- suppressWarnings(as.numeric(pgs$PGS))

ph <- read.table("pheno_audpc.txt", header=FALSE, stringsAsFactors=FALSE)
if (nrow(ph) > 0 && (grepl("FID", ph[1,1], ignore.case=TRUE) || grepl("IID", ph[1,2], ignore.case=TRUE))) {
  ph <- ph[-1, , drop=FALSE]
}
ph <- ph[,1:3]
names(ph) <- c("FID","IID","AUDPC")
ph$AUDPC <- suppressWarnings(as.numeric(ph$AUDPC))

m <- merge(pgs, ph, by=c("FID","IID"), all.x=TRUE)
write.table(m, "gs-lite/pgs_with_audpc.tsv", sep="\t", row.names=FALSE, quote=FALSE)
mr <- m[order(-m$PGS), ]
write.table(mr, "gs-lite/pgs_ranked.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(head(mr, 20), "gs-lite/top20_by_PGS.tsv", sep="\t", row.names=FALSE, quote=FALSE)

png("gs-lite/PGS_vs_AUDPC.png", width=1000, height=800, res=150)
plot(m$PGS, m$AUDPC, xlab="Polygenic Score (PGS)", ylab="AUDPC",
     main="PGS vs AUDPC (quick check)")
abline(lm(AUDPC ~ PGS, data=m), col="black")
r <- suppressWarnings(cor(m$PGS, m$AUDPC, use="complete.obs"))
if (!is.na(r)) legend("topleft", bty="n", legend=paste0("r = ", round(r, 3)))
dev.off()
EOF

echo "[GS-lite] Done."
echo "[GS-lite] Results in: $OUTDIR/"
```

Submit the job:
```bash
sbatch 50_gs_lite_audpc.sh
```


### Step 2 — Inspect outputs

```bash
# The per-line PGS file (look for a SCORE* column)
head -n 30 gs-lite/gs_pgs.profile

# Ranked table of lines by PGS (descending)
head -n 30 gs-lite/pgs_ranked.tsv

# Top 20 lines by PGS
cat gs-lite/top20_by_PGS.tsv

# Quick check of the plot
ls -lh gs-lite/PGS_vs_AUDPC.png
```

**Files produced**  

- `gs-lite/gs_pgs.profile` — PLINK’s per-individual score file (has one `SCORE*` column).
- `gs-lite/pgs_with_audpc.tsv` — merged table with `FID IID PGS AUDPC` (unsorted).
- `gs-lite/pgs_ranked.tsv` — same, **sorted by PGS** (highest first).
- `gs-lite/top20_by_PGS.tsv` — top 20 individuals by PGS.
- `gs-lite/PGS_vs_AUDPC.png` — scatter + regression line, with correlation `r` (if computable).
- `gs-lite/gs_pgs.log` — PLINK log for the scoring step.
- If your BIM needed position-based IDs: `gs-lite/data_posids.*` — an auto-converted copy used by the script.


### Step 3 - Interpreting results for breeding decisions

- **PGS sign and magnitude:** A **higher PGS** indicates the individual carries more alleles with **positive BETA** for **AUDPC**. Depending on how AUDPC is defined in your experiment (larger values = **more disease**), you’ll either select **low-PGS** (if high AUDPC is bad) or **high-PGS** (if high AUDPC is good).  
  - *Most disease contexts:* **Lower AUDPC is better** → **prefer lower PGS** (if BETA was fit to raw AUDPC).  
  - If you transformed AUDPC (e.g., inverse-normal) or reversed the trait, align interpretation accordingly.
- **Ranking vs. phenotype:** If the **PGS vs. AUDPC plot** shows a clear relationship (|r| appreciable), even a small marker set gives you **triage power**: shortlist candidates **before** full phenotyping.
- **Genomic selection mindset:** This GS-lite uses only a handful of GWAS SNPs. For **production GS**, you’ll typically use **many thousands** of markers and a **predictive model** (e.g., GBLUP/RR-BLUP/Bayesian) with **train/test** partitions or cross-validation.
- **Stacking with other traits:** You can compute separate PGS for **yield/quality** traits, then make a **multi-trait index** (weighted sum) to balance **disease** vs **productivity**.
- **Environment and G×E:** If disease pressure varies with **heat** or **humidity**, you may re-estimate BETA in environment-specific subsets or include **weather covariates**, then recompute PGS for environment-tailored selection.

**Actionable summary:**  
1. Use `pgs_ranked.tsv` to identify the **top-N candidates** consistent with your breeding goal (e.g., **lowest PGS** if lower AUDPC is better).  
2. Cross-reference with **field performance** and **other key traits**.  
3. Consider **validation** in an independent block/year.  
4. Graduate to a **full GS model** as soon as you have enough training data.

#### Notes / Gotchas  
- If **no Bonferroni hits**, the script falls back to **top-20 p-values** (or top-50 from the GWAS file). This still yields a teaching-ready PGS but expect a weaker correlation.
- If your `.bim` has many `.` SNP IDs, the script **auto-relabels** them to `CHR:BP` and rewrites the weights to match, ensuring PLINK’s `--score` can find them.
- Always confirm that the **effect allele (A1)** in the weights matches the allele coding in your GWAS file (we use your tables’ `A1` directly).


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
