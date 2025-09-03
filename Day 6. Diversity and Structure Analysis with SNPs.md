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

The dataset we will use comes from the study by Hazzouri *et al.* (2019), published in *Nature Communications*:  
*Genome-wide association mapping of date palm fruit traits.*  
**Link:** (https://www.nature.com/articles/s41467-019-12604-9)

In this work, the authors:  
- Produced a **high-quality genome assembly** for the date palm (*Phoenix dactylifera*).  
- Conducted **GWAS** on fruit traits, identifying candidate genes such as the **R2R3-MYB transcription factor VIRESCENS** (fruit color) and **invertases** (sugar composition).  
- Generated a large SNP dataset across hundreds of individuals.  

For training purposes, I prepared a **reduced version** of their SNP dataset:  

- **File**: `dataset120_chr18.vcf.gz`  
- **Individuals**: Subset of **120 individuals** from the original population.  
- **Variants**: Includes SNPs from **contigs 1–18 only**, with a **reduced number of SNPs**.  

This reduced VCF is representative and enables us to run the full pipeline — including **PCA, ADMIXTURE, LD pruning/decay, and diversity statistics** — in a practical timeframe.


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
   #SBATCH --mem=1G
   #SBATCH --time=00:30:00
   #SBATCH -o vcf2plink.out
   #SBATCH -e vcf2plink.err

   module load PLINK

   # Input VCF (from Day 5)
   IN_VCF="/lisc/data/scratch/course/pgbiow/data/VCF/dataset120_chr18.vcf.gz"

   echo "Converting ${IN_VCF} to PLINK binary format..."
   plink --vcf "${IN_VCF}" \
         --make-bed \
         --allow-extra-chr \
         --double-id \
         --out ./plink/my_data

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
   #SBATCH --mem=1G
   #SBATCH --time=00:30:00
   #SBATCH -o ld_prune.out
   #SBATCH -e ld_prune.err

   module load PLINK

   IN_BASE="./plink/my_data"  # from Step 0.1

   echo "Selecting approximately independent SNPs (LD pruning)..."
   plink --bfile "${IN_BASE}" \
         --allow-extra-chr \
         --indep-pairwise 50 5 0.2 \
         --out ./plink/my_data_prune

   echo "Creating pruned dataset..."
   plink --bfile "${IN_BASE}" \
         --allow-extra-chr \
         --extract ./plink/my_data_prune.prune.in \
         --make-bed \
         --out ./plink/my_data_pruned

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

Open the `.log` file to check how many SNPs were removed after filtering:
   ```bash
   cat ld_prune.out
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

   IN_BASE="./plink/my_data"  # unpruned dataset from Step 0.1

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

   IN_VCF=/lisc/data/scratch/course/pgbiow/data/VCF/dataset120_chr18.vcf.gz"  # same input used in Step 0.1

   echo "Per-individual heterozygosity and inbreeding coefficient..."
   vcftools --gzvcf "${IN_VCF}" \
            --het \
            --out ./diversity/vcftools_het

   echo "Per-individual missingness..."
   vcftools --gzvcf "${IN_VCF}" \
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

### Step 3 — FST Between Two Populations (VCFtools)

We can calculate **Weir & Cockerham’s FST** between two populations using VCFtools.


**1. Create population files**  

From the metadata CSV, extract the sample IDs for each population. Each file must contain **one sample ID per line**.

```bash
# Metadata file
META=/lisc/data/scratch/course/pgbiow/data/metadata/gwas_pop_table_120.csv

# Create pop1.txt (example: Al-Ain - Abu Dhabi)
awk -F',' 'NR>1 {gsub(/"/,""); if($2=="Al-Ain - Abu Dhabi") print $1}' $META > ./diversity/pop1.txt

# Create pop2.txt (replace with the second population name)
awk -F',' 'NR>1 {gsub(/"/,""); if($2=="Ras al-Khaimah") print $1}' $META > ./diversity/pop2.txt

# Check that both lists have samples
wc -l ./diversity/pop1.txt ./diversity/pop2.txt
```

**2. Create a script 12_vcftools_fst.sh:**

1. Open:
   ```bash
   vi 12_vcftools_fst.sh
   ```
2. Press `i`.  
3. **Copy & paste**:
   
```bash
#!/bin/bash
#SBATCH --job-name=vcftools_fst
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH -o diversity/vcftools_fst.out
#SBATCH -e diversity/vcftools_fst.err

module load vcftools

vcftools --gzvcf /lisc/data/scratch/course/pgbiow/data/VCF/dataset120_chr18.vcf.gz \
         --weir-fst-pop ./diversity/pop1.txt \
         --weir-fst-pop ./diversity/pop2.txt \
         --out ./diversity/fst_result
```
4. To save: `Esc`, `:wq`, `Enter`.  

5. Submit:
```bash
sbatch 12_vcftools_fst.sh
```
6. Inspect results:
```bash
head ./diversity/fst_result.weir.fst
```

To compute the mean FST:
```bash
awk 'NR>1 && $3!="nan"{s+=$3; n++} END{if(n>0) print "Mean FST = " s/n; else print "No valid sites"}' ./diversity/fst_result.weir.fst
```

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
   INPUT_BASE_PRUNED="./plink/my_data_pruned"

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


```r
# ================================
# PCA Visualization Script
# ================================

# Load library
library(ggplot2)

# --- Load PCA results from PLINK ---
eigenvec <- read.table("pca_results.eigenvec", header = FALSE)
colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:(ncol(eigenvec)-2)))

# Load eigenvalues
eigenval <- scan("pca_results.eigenval")

# Calculate percentage variance explained
variance_explained <- round(100 * eigenval / sum(eigenval), 2)

# Print the variance explained for the first 5 PCs
cat("Variance explained (%):\n")
print(variance_explained[1:5])

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

### Step 1 - Run ADMIXTURE for different K values

ADMIXTURE takes PLINK BED files as input. We will run it for `K=2` to `K=5` as an example. It is recommended to run ADMIXTURE multiple times for each K value with different random seeds to check for convergence.

**Create a file named `run_admixture.sh` using `vi`:**

```bash
#!/bin/bash
#SBATCH --job-name=run_admixture
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH -o admixture.out
#SBATCH -e admixture.err

module load admixture

# Input dataset: LD-pruned files created in Step 2 (my_data_pruned.bed/.bim/.fam)
INPUT_BASE_PRUNED="my_data_pruned"

# Loop through different K values
for K in {2..5}
do
    echo "Running ADMIXTURE for K=$K"
    admixture --cv -j4 $INPUT_BASE_PRUNED.bed $K | tee admixture_K${K}.log
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

### Step 2 - Determine Optimal K and Visualize Admixture Results (Conceptual)

To find the optimal `K`, you typically look for the `K` value with the lowest cross-validation error (CV error) from the `.log` files generated by ADMIXTURE.

```bash
# Extract CV errors from log files
grep "CV error" admixture_K*.log
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
# ================================
# ADMIXTURE Visualization Script
# ================================

# Load library
library(ggplot2)

# --- Step 1: Load Q file (change K value as needed) ---
# Example for K=3
Qfile <- "my_data_pruned.3.Q"

q_df <- read.table(Qfile, header = FALSE)
colnames(q_df) <- paste0("Cluster", 1:ncol(q_df))

# --- Step 2: Add individual IDs ---
# Read sample IDs from PLINK .fam file
fam <- read.table("my_data_pruned.fam")
q_df$IID <- fam$V2   # Second column is individual ID

# --- Step 3: Reshape data for ggplot ---
library(reshape2)
q_melt <- melt(q_df, id.vars = "IID", variable.name = "Cluster", value.name = "Ancestry")

# --- Step 4: Plot barplot ---
ggplot(q_melt, aes(x = IID, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
  labs(title = "ADMIXTURE Results (K=3)",
       x = "Individuals",
       y = "Ancestry Proportion") +
  scale_fill_brewer(palette = "Set2")
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


### Step 2 — Visualize SNP Density Plot in R

Now that we have `snp_density.txt`, we can plot SNP density across the genome using **R**.

#### Example R Script

Copy this into your R session (on the cluster with graphics enabled or in RStudio):

```r
# ================================
# SNP Density Plot Script
# ================================

# Load libraries
library(ggplot2)
library(dplyr)

# --- Step 1: Load SNP density data ---
# snp_density.txt has 4 columns:
# count  chromosome  window_start  window_end
snp_df <- read.table("snp_density.txt", header = FALSE)
colnames(snp_df) <- c("Count", "Chromosome", "Window_Start", "Window_End")

# --- Step 2: Order chromosomes (if numeric, convert properly) ---
snp_df$Chromosome <- factor(snp_df$Chromosome,
                            levels = mixedsort::mixedsort(unique(snp_df$Chromosome)))

# --- Step 3: Make plot ---
ggplot(snp_df, aes(x = Window_Start, y = Count)) +
  geom_line(color = "steelblue") +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  theme_minimal(base_size = 14) +
  labs(title = "SNP Density Across the Genome",
       x = "Genomic Position (bp)",
       y = "Number of SNPs per Window") +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
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
