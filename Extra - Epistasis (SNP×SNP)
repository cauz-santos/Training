# Epistasis (SNP×SNP) for SUC

**Goal:** Test whether pairs of GWAS-relevant SNPs interact (depart from additivity) for the quantitative trait **SUC**, controlling for **PC1–PC5**.  

---

## 6.1 — Pick candidate SNPs (limit the search space)

Take the **top 300 additive GWAS SNPs by p-value** (change `N_TOP` as you like) to keep pair counts manageable.

```bash
vi 60_pick_epistasis_candidates.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=epi_pick
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:02:00
#SBATCH -o epi_pick.out
#SBATCH -e epi_pick.err
set -euo pipefail

IN="./gwas/gwas_suc_linear.assoc.linear"   # from Part 2
OUTDIR="epistasis"; mkdir -p "$OUTDIR"
N_TOP=300

# Robust header-aware extraction of top N additive SNPs by p-value
awk 'BEGIN{FS=OFS="\t"}
NR==1{
  for(i=1;i<=NF;i++){h[$i]=i}
  # Require these columns to be present
  if(!("TEST" in h) || !("P" in h) || !("SNP" in h)){
    print "ERROR: Required columns TEST/P/SNP not found" > "/dev/stderr"; exit 1
  }
  next
}
$h["TEST"]=="ADD" && $h["P"]!="" && $h["P"]!="NA" {
  print $h["SNP"], $h["P"]
}' "$IN" | sort -g -k2,2 | head -n ${N_TOP} | cut -f1 > "$OUTDIR/epi_candidates.snps"

echo "Wrote: $OUTDIR/epi_candidates.snps (top ${N_TOP} SNPs by p-value)"
```

```bash
sbatch 60_pick_epistasis_candidates.sh
```

> **Tip:** With ~120 samples, your QC already used `--maf 0.05`, which helps avoid tiny genotype cells in the 3×3 grid per SNP (0/1/2).

---

## 6.2 — Run PLINK epistasis (linear model + PCs)

PLINK’s `--epistasis` fits a linear model for SUC **including the main effects** and the **SNP×SNP interaction**, while PCs enter as covariates. We also set a print gate (`--epi1 1e-4`) so files don’t explode.

```bash
vi 61_run_epistasis_linear.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=epi_run
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH -o epi_run.out
#SBATCH -e epi_run.err

module load PLINK

BFILE="./plink/gwas_data_qc"
PHENO="pheno_suc.txt"
COVAR="covar_pcs.txt"
SNPLIST="epistasis/epi_candidates.snps"

plink --bfile "$BFILE" \
      --extract "$SNPLIST" \
      --pheno "$PHENO" --pheno-name SUC \
      --covar "$COVAR" --covar-name PC1,PC2,PC3,PC4,PC5 \
      --allow-extra-chr --allow-no-sex \
      --epistasis \
      --epi1 1e-4 \
      --out epistasis/epi_suc_linear

echo "Done: epistasis/epi_suc_linear.epi.qt"
```

```bash
sbatch 61_run_epistasis_linear.sh
```

**Output:** `epistasis/epi_suc_linear.epi.qt`  
(Contains `SNP1 SNP2 ... P`, where **P** is the interaction p-value for the SNP×SNP term.)

---

## 6.3 — Summarize + FDR and make a 9-cell plot for the top pair

This R script:
1) Reads the epistasis results, computes **FDR (BH)**, writes `top50_epistasis.tsv`.  
2) Extracts genotypes (dosage 0/1/2) for the **top pair** via PLINK `--recode A`.  
3) Plots the **9 genotype-combination cells** with mean SUC and sample counts.

```bash
vi 62_epistasis_postprocess.R
```

```r
# Epistasis post-processing + 9-cell plot
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# 1) Read epistasis table and compute FDR
epi <- fread("epistasis/epi_suc_linear.epi.qt")
# Standardize colnames and check essentials
names(epi) <- toupper(names(epi))
stopifnot(all(c("SNP1","SNP2","P") %in% names(epi)))
epi <- epi[!is.na(P)]
epi[, FDR := p.adjust(P, "BH")]
setorder(epi, P)
fwrite(head(epi, 50), "epistasis/top50_epistasis.tsv", sep = "\t")

lead <- epi[1]
cat(sprintf("Top pair: %s × %s | P=%.3e | FDR=%.3e\n", lead$SNP1, lead$SNP2, lead$P, lead$FDR))

# 2) Get genotypes (0/1/2) for the top pair
cmd <- sprintf(
  'plink --bfile plink/gwas_data_qc --snps %s,%s --recode A --out epistasis/top_pair --allow-extra-chr --allow-no-sex',
  lead$SNP1, lead$SNP2
)
system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)

geno <- fread("epistasis/top_pair.raw")   # columns: FID IID and SNP dosage columns
# Locate exact genotype columns (handle possible name quirks)
cn <- names(geno)
g1col <- cn[cn == lead$SNP1]
g2col <- cn[cn == lead$SNP2]
# Fallback: prefix match if exact not found
if (!length(g1col)) g1col <- cn[grep(paste0("^", gsub("([][(){}.+*?^$|\\\\])","\\\\\\1", lead$SNP1)), cn)]
if (!length(g2col)) g2col <- cn[grep(paste0("^", gsub("([][(){}.+*?^$|\\\\])","\\\\\\1", lead$SNP2)), cn)]
stopifnot(length(g1col)>=1, length(g2col)>=1)
g1col <- g1col[1]; g2col <- g2col[1]

ph <- fread("pheno_suc.txt", header = FALSE)
setnames(ph, c("FID","IID","SUC"))

m <- merge(
  geno[, .(FID, IID, g1 = get(g1col), g2 = get(g2col))],
  ph, by = c("FID","IID"), all.x = TRUE
)

# 3) 9-cell grid (g1 x g2), with mean SUC and counts
m[, g1 := as.integer(round(g1))]
m[, g2 := as.integer(round(g2))]
cell <- m[!is.na(SUC), .(mean_SUC = mean(SUC), n = .N), by = .(g1, g2)]

# Warn if any cell is very small (e.g., n<5)
small_cells <- cell[n < 5]
if (nrow(small_cells) > 0) {
  msg <- paste0("Warning: small cells detected (n<5) in genotype grid: ",
                paste(apply(small_cells, 1, function(r) sprintf("(%s,%s):n=%s", r["g1"], r["g2"], r["n"])), collapse = "; "))
  message(msg)
}

# Heatmap of mean SUC with counts
p <- ggplot(cell, aes(x = factor(g1), y = factor(g2), fill = mean_SUC,
                      label = paste0(round(mean_SUC, 2), "\\n(n=", n, ")"))) +
  geom_tile() +
  geom_text(size = 3) +
  scale_fill_continuous(name = "Mean SUC") +
  labs(x = paste0(lead$SNP1, " genotype (0/1/2)"),
       y = paste0(lead$SNP2, " genotype (0/1/2)"),
       title = sprintf("Epistasis 9-cell plot: %s × %s (P=%.2e, FDR=%.2e)",
                       lead$SNP1, lead$SNP2, lead$P, lead$FDR)) +
  theme_minimal(base_size = 12)

ggsave("epistasis/epistasis_9cell_toppair.png", p, width = 6, height = 5, dpi = 150)
cat("Saved: epistasis/epistasis_9cell_toppair.png\n")
```

```bash
Rscript 62_epistasis_postprocess.R
```

**Outputs created:**  
- `epistasis/top50_epistasis.tsv` — top pairs with raw P and FDR.  
- `epistasis/epistasis_9cell_toppair.png` — heatmap of mean SUC across genotype combinations.

---

- **Interpretation:** The 9-cell heatmap; if the pattern deviates from a plane (e.g., curvature), that supports interaction.  
- **QC sanity:** ensure both SNPs have **MAF ≥ 0.05** and avoid tiny cells (e.g., any with `n<5`).  
- **Multiple testing:** you limited the search (top-SNP set) and applied **FDR**; for publication-grade evidence, **replicate** in an independent set.  
- **Breeding angle:** interacting loci can yield **context-dependent effects**; consider haplotypes or two-marker assays during validation.

