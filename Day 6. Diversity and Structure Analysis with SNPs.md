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

We will work with a VCF file already filtered (e.g., for depth, quality, missingness, and minor allele frequency).

- `my_filtered_variants.vcf.gz` — High-quality, biallelic SNP dataset

This file will be converted into formats required for **PLINK**, **ADMIXTURE**, and other tools.



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

We will convert a **filtered VCF** into PLINK format and perform **LD pruning**. The pruned dataset will be used for **PCA** (and later for ADMIXTURE/GWAS).

### Step 1 — Convert VCF to PLINK format

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
   #SBATCH --job-name=vcf2plink
   #SBATCH --cpus-per-task=1
   #SBATCH --mem=8G
   #SBATCH --time=00:30:00
   #SBATCH -o vcf2plink.out
   #SBATCH -e vcf2plink.err

   module load plink

   # Input VCF (from Day 5)
   IN_VCF="my_filtered_variants.vcf.gz"

   echo "Converting ${IN_VCF} to PLINK binary format..."
   plink --vcf "${IN_VCF}" \
         --make-bed \
         --out my_data

   echo "Done. Created: my_data.bed / .bim / .fam"
   ```
4. Press `Esc`, type `:wq` and press `Enter` to save and exit.  
5. Submit the job:
   ```bash
   sbatch 00_convert_vcf_to_plink.sh
   ```


### Step 2 — LD Pruning (remove linked SNPs)

**What is LD Pruning?**  

**Linkage Disequilibrium (LD)** refers to the non-random association of alleles at different loci.  
In simple terms, if two SNPs are very close to each other on the chromosome, they are often inherited together — meaning their genotypes are **correlated**.  

- Example: If SNP A and SNP B are always observed together in your dataset, they are in **high LD**.  
- For analyses like **PCA** or **Admixture**, including both SNPs does not add new information — it only adds redundancy.  

**Why do we prune SNPs in LD?**  

- **PCA & Admixture assume independence**: correlated SNPs inflate the signal and may bias the results.  
- **Computational efficiency**: fewer SNPs make analyses faster without losing meaningful information.  
- **Interpretability**: a pruned dataset captures the true population structure instead of local chromosomal effects.  

**How does pruning work in PLINK?**  

PLINK uses a **sliding window** approach:
- It scans windows of a defined number of SNPs (e.g., 50 SNPs).  
- Within each window, it calculates the **pairwise correlation (r²)** between SNPs.  
- If two SNPs are too correlated (above a threshold, e.g., r² ≥ 0.2), one of them is removed.  
- The window slides forward (e.g., by 5 SNPs), and the process repeats until the genome is scanned.  


**Create the Slurm script with `vi`:**

1. Open the editor:
   ```bash
   vi 01_ld_pruning.sh
   ```
2. Press `i` to enter *INSERT* mode.  
3. **Copy & paste**:
   ```bash
   #!/bin/bash
   #SBATCH --job-name=ld_prune
   #SBATCH --cpus-per-task=1
   #SBATCH --mem=8G
   #SBATCH --time=00:30:00
   #SBATCH -o ld_prune.out
   #SBATCH -e ld_prune.err

   module load plink

   IN_BASE="my_data"  # from Step 0.1

   echo "Selecting approximately independent SNPs (LD pruning)..."
   plink --bfile "${IN_BASE}" \
         --indep-pairwise 50 5 0.2 \
         --out my_data_prune

   echo "Creating pruned dataset..."
   plink --bfile "${IN_BASE}" \
         --extract my_data_prune.prune.in \
         --make-bed \
         --out my_data_pruned

   echo "Done. Use 'my_data_pruned' for PCA and structure analyses."
   ```
4. Press `Esc`, then `:wq`, `Enter`.  
5. Submit the job:
   ```bash
   sbatch 01_ld_pruning.sh
   ```

**Outputs created**

- `my_data_prune.prune.in` — SNPs to keep  
- `my_data_prune.prune.out` — SNPs removed  
- `my_data_pruned.*` — **LD-pruned** PLINK dataset (use this for **PCA**)

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
   #SBATCH --mem=8G
   #SBATCH --time=00:20:00
   #SBATCH -o plink_div.out
   #SBATCH -e plink_div.err

   module load plink

   IN_BASE="my_data"  # unpruned dataset from Step 0.1

   echo "Calculating heterozygosity and inbreeding coefficient (per sample)..."
   plink --bfile "${IN_BASE}" \
         --het \
         --out plink_het

   echo "Calculating per-individual missingness..."
   plink --bfile "${IN_BASE}" \
         --missing \
         --out plink_missing

   echo "Done. Outputs: plink_het.het, plink_missing.imiss"
   ```
4. `Esc`, `:wq`, `Enter`.  
5. Submit:
   ```bash
   sbatch 10_plink_diversity.sh
   ```

**How to read the outputs**

- `plink_het.het` → O(HOM), E(HOM), N_SNPs, **F** = (E - O) / E  
- `plink_missing.imiss` → fraction of missing genotypes per individual


### Step 2 — VCFtools: MAF, Heterozygosity, Missingness  

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
   #SBATCH --mem=6G
   #SBATCH --time=00:20:00
   #SBATCH -o vcftools_div.out
   #SBATCH -e vcftools_div.err

   module load vcftools

   IN_VCF="my_filtered_variants.vcf.gz"  # same input used in Step 0.1

   echo "Per-site allele frequencies (MAF)..."
   vcftools --gzvcf "${IN_VCF}" \
            --freq \
            --out vcftools_maf

   echo "Per-individual heterozygosity and inbreeding coefficient..."
   vcftools --gzvcf "${IN_VCF}" \
            --het \
            --out vcftools_het

   echo "Per-individual missingness..."
   vcftools --gzvcf "${IN_VCF}" \
            --missing-indv \
            --out vcftools_missing

   echo "Done. Outputs: vcftools_maf.frq, vcftools_het.het, vcftools_missing.imiss"
   ```
4. `Esc`, `:wq`, `Enter`.  
5. Submit:
   ```bash
   sbatch 11_vcftools_diversity.sh
   ```

**Key outputs**

- `vcftools_maf.frq` → allele frequencies per SNP (inspect MAF distribution)  
- `vcftools_het.het` → per-sample O(HOM), E(HOM), **F**  
- `vcftools_missing.imiss` → per-sample missingness


### Diversity QC — Questions to Consider

- Are there individuals with **very high missingness** (e.g., >10–20%)?  
- Do some samples show **unusually low heterozygosity** (possible inbreeding or data issues)?  
- Is the **MAF distribution** dominated by very rare variants (might affect power/interpretation)?

---

## Part 3 — Principal Component Analysis (PCA) with PLINK
Principal Component Analysis (PCA) is a statistical method that reduces the dimensionality of complex datasets while retaining most of the variation. In population genetics, PCA is used to visualize genetic relationships between individuals and identify population structure. Individuals that are genetically similar will cluster together in the PCA plot, while genetically distinct populations will separate.

We will run PCA on the **LD-pruned** dataset produced in **Step 2** (`my_data_pruned.*`).  

**Inputs used here**

- `my_data_pruned.bed/.bim/.fam` (created by `01_ld_pruning.sh`)

**Outputs produced**

- `pca_results.eigenval` — variance explained by each PC  
- `pca_results.eigenvec` — PC scores per individual  

### Step 2.1 — Run PCA (create script with `vi`)

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
   #SBATCH --mem=8G
   #SBATCH --time=00:30:00
   #SBATCH -o run_pca.out
   #SBATCH -e run_pca.err

   module load plink

   # LD-pruned dataset from Part 0
   INPUT_BASE_PRUNED="my_data_pruned"

   echo "Running PCA on ${INPUT_BASE_PRUNED}..."
   plink --bfile "${INPUT_BASE_PRUNED}" \
         --pca \
         --out pca_results

   echo "Done. Outputs: pca_results.eigenval, pca_results.eigenvec"
   ```
4. `Esc`, `:wq`, `Enter`.  
5. Submit:
   ```bash
   sbatch 20_run_pca.sh
   ```

## 📊 PCA Visualization in R

After running PCA with PLINK, we will visualize the results using R.  
Copy the following script into RStudio (or an interactive R session on the cluster) and run it.

```r
# ================================
# PCA Visualization Script
# ================================

# Load library
library(ggplot2)

# --- Load PCA results from PLINK ---
eigenvec <- read.table("pca_results.eigenvec", header = FALSE)
colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:(ncol(eigenvec)-2)))

eigenval <- scan("pca_results.eigenval")
variance_explained <- round(100 * eigenval / sum(eigenval), 2)

# --- Load population labels ---
# File must contain at least two columns: IID and Population
# Example:
# IID    Population
# sample1   PopA
# sample2   PopB
popinfo <- read.table("population_labels.txt", header = TRUE, stringsAsFactors = FALSE)

# --- Merge data ---
pca_df <- merge(eigenvec, popinfo, by = "IID")

# --- Plot PCA ---
ggplot(pca_df, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA of Genomic Data (colored by population)",
       x = paste0("PC1 (", variance_explained[1], "% variance)"),
       y = paste0("PC2 (", variance_explained[2], "% variance)")) +
  scale_color_brewer(palette = "Set1")
```

> This script will open a PCA scatterplot window with individuals colored by population.


**Interpreting PCA**

- **PC1, PC2** typically separate major population groups; nearby individuals in the plot are genetically similar.  
- Use `pca_results.eigenval` to compute the **variance explained** (e.g., PC1% = eigenval₁ / sum(eigenvals) × 100).  
- The file `pca_results.eigenvec` begins with `FID IID` followed by PC1, PC2, … per individual.

---

## Part 4: Admixture Analysis

**What is Admixture?**

ADMIXTURE is a software tool for estimating individual ancestries from multilocus SNP genotype datasets. It assumes a model where individuals originate from a mixture of `K` ancestral populations. The output provides the proportion of ancestry each individual derives from each of the `K` ancestral populations. Unlike PCA, which is an unsupervised method, Admixture provides a more direct estimate of ancestral contributions.

We will run ADMIXTURE for different values of `K` (number of ancestral populations) and then determine the optimal `K`.

### 3.1: Run ADMIXTURE for different K values

ADMIXTURE takes PLINK BED files as input. We will run it for `K=2` to `K=5` as an example. It is recommended to run ADMIXTURE multiple times for each K value with different random seeds to check for convergence.

**Create a file named `run_admixture.sh` using `vi`:**

```bash
#!/bin/bash
#SBATCH --job-name=run_admixture
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH -o admixture.out
#SBATCH -e admixture.err

# Load necessary modules
module load admixture

# Define input PLINK files base name (pruned data)
INPUT_BASE_PRUNED="my_data_pruned"

# Loop through different K values
for K in {2..5}
do
    echo "Running ADMIXTURE for K=$K"
    admixture --cv $INPUT_BASE_PRUNED.bed $K | tee admixture_K${K}.log
done

echo "ADMIXTURE runs complete."
```

**Submit the job:**

```bash
sbatch run_admixture.sh
```

*   The `--cv` flag performs cross-validation, which helps in determining the optimal `K` value.
*   The `tee` command saves the output to both the screen and a log file.

**❓ Question:** What does the `K` parameter in ADMIXTURE represent? How does the algorithm estimate ancestry proportions?

### 3.2: Determine Optimal K and Visualize Admixture Results (Conceptual)

To find the optimal `K`, you typically look for the `K` value with the lowest cross-validation error (CV error) from the `.log` files generated by ADMIXTURE.

```bash
# Extract CV errors from log files
grep "CV error" admixture_K*.log
```

**Conceptual Visualization:**

Admixture results are often visualized as stacked bar plots, where each bar represents an individual, and segments within the bar represent the proportion of ancestry from each ancestral population. This visualization is usually done in R or Python.

```python
# Conceptual Python code for Admixture visualization
# You would need to transfer the .Q files (e.g., my_data_pruned.2.Q) to your local machine
# and use a plotting environment.

# import pandas as pd
# import matplotlib.pyplot as plt

# Load Q file (e.g., for K=2)
# q_df = pd.read_csv("my_data_pruned.2.Q", sep="\s+", header=None)

# Create a dummy individual ID column if not present
# q_df["Individual"] = [f"Ind{i+1}" for i in range(len(q_df))]

# Plotting
# q_df.set_index("Individual").plot(kind="bar", stacked=True, figsize=(12, 6), legend=False)
# plt.title("Admixture Proportions (K=2)")
# plt.xlabel("Individual")
# plt.ylabel("Ancestry Proportion")
# plt.xticks(rotation=90, fontsize=8)
# plt.tight_layout()
# plt.show()
```

**❓ Question:** How do you interpret a stacked bar plot from ADMIXTURE? What does it mean if an individual has ancestry from multiple `K` populations?

---

## Part 4: SNP Density Plots

**What are SNP Density Plots?**

SNP density plots visualize the distribution of Single Nucleotide Polymorphisms (SNPs) across the genome. They help identify regions with high or low genetic variation, which can be indicative of selective sweeps, recombination hotspots, or genomic rearrangements. These plots are typically generated by counting SNPs within fixed-size genomic windows.

### 4.1: Prepare Data for SNP Density Plot

We need to extract chromosome and position information from our VCF file and then bin the SNPs into genomic windows. We can use `bcftools query` and some shell scripting for this.

**Create a file named `snp_density_data.sh` using `vi`:**

```bash
#!/bin/bash
#SBATCH --job-name=snp_density_data
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH -o snp_density_data.out
#SBATCH -e snp_density_data.err

# Load necessary modules
module load bcftools

# Define input VCF file (use the filtered VCF)
INPUT_VCF="my_filtered_variants.vcf.gz"

# Define window size (e.g., 100 kb = 100000 bp)
WINDOW_SIZE=100000

echo "Extracting SNP positions and calculating density..."

# Extract CHROM and POS, then calculate window for each SNP
bcftools query -f 
'[%CHROM\t%POS\n]
' $INPUT_VCF | \
awk -v WS=$WINDOW_SIZE 
'BEGIN {OFS="\t"} {print $1, int($2/WS)*WS, int($2/WS)*WS + WS - 1}
' | \
sort -k1,1 -k2,2n | \
uniq -c > snp_density.txt

echo "SNP density data generated: snp_density.txt"
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

**❓ Question:** How does changing the `WINDOW_SIZE` affect the appearance and interpretation of a SNP density plot?

### 4.2: Visualize SNP Density Plot (Conceptual)

SNP density plots are typically line plots or bar plots showing the number of SNPs per window across chromosomes. This is best done in R or Python.

```python
# Conceptual Python code for SNP Density Plot visualization
# You would need to transfer snp_density.txt to your local machine
# and use a plotting environment.

# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns

# Load the SNP density data
# snp_density_df = pd.read_csv("snp_density.txt", sep="\s+", header=None, names=["Count", "CHROM", "Window_Start", "Window_End"])

# Sort by chromosome and window start for proper plotting order
# snp_density_df["CHROM_numeric"] = snp_density_df["CHROM"].astype("category").cat.codes # For consistent sorting
# snp_density_df = snp_density_df.sort_values(by=["CHROM_numeric", "Window_Start"])

# Create a combined genomic position for plotting across chromosomes
# snp_density_df["Genomic_Position"] = snp_density_df.groupby("CHROM")["Window_Start"].transform(lambda x: x - x.min())
# cumulative_lengths = snp_density_df.groupby("CHROM")["Window_End"].max().cumsum().shift(1).fillna(0)
# snp_density_df["Genomic_Position"] = snp_density_df.apply(lambda row: row["Genomic_Position"] + cumulative_lengths[row["CHROM"]], axis=1)

# Plotting
# plt.figure(figsize=(15, 6))
# sns.lineplot(data=snp_density_df, x="Genomic_Position", y="Count", hue="CHROM", palette="tab10", linewidth=0.8)
# plt.title("SNP Density Across the Genome (100kb Windows)")
# plt.xlabel("Genomic Position")
# plt.ylabel("Number of SNPs")
# plt.xticks([]) # Hide x-axis ticks for clarity if many chromosomes
# plt.legend(title="Chromosome", bbox_to_anchor=(1.05, 1), loc="upper left")
# plt.tight_layout()
# plt.show()
```

**❓ Question:** What patterns in a SNP density plot might suggest regions under positive selection or regions with low recombination rates?

---

## Conclusion

Congratulations! You have completed a comprehensive practical on population structure analysis and SNP density plotting. You have learned how to prepare VCF data for specialized software, perform PCA and Admixture analysis to infer population relationships and ancestries, and generate SNP density plots to visualize genomic variation. These are powerful techniques for understanding the evolutionary history and genetic diversity of populations.

Remember that these analyses are often iterative, and the interpretation of results requires careful consideration of biological context and potential biases. Keep exploring and applying these tools to your own research!

