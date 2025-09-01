## Week 2: Applied Bioinformatics for Genomics and Breeding

## Day 6: Diversity and Structure Analysis with SNPs

Welcome to **Day 6** of our bioinformatics training!  

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

We will work with a VCF file already filtered during Day 5 (e.g., for depth, quality, missingness, and minor allele frequency).

- `my_filtered_variants.vcf.gz` — High-quality, biallelic SNP dataset

This file will be converted into formats required for **PLINK**, **ADMIXTURE**, and other tools.



### Why Does This Matter?

Understanding **genetic diversity and structure** is essential in genomics and breeding:

- In **breeding programs**, it helps identify divergent individuals, avoid inbreeding, and guide parent selection.
- In **population genetics**, it reveals how populations are related and how they evolved.
- In **conservation**, it supports decisions on which populations or individuals to prioritize.

The tools we explore today — like PCA and ADMIXTURE — are **foundational methods** in modern genomics. Combined with diversity metrics (e.g., heterozygosity, allele frequency spectra), they allow you to **quantify and visualize genetic patterns** in your dataset.

---

**Table 1. Tools You Will Use**

| Tool       | Purpose                                      |
|------------|----------------------------------------------|
| **PLINK**      | Data conversion, PCA, diversity stats        |
| **ADMIXTURE**  | Ancestry inference (admixture proportions)  |
| **VCFtools**   | Summary statistics, missingness, MAF        |
| **BCFtools**   | SNP extraction, VCF querying                 |
| **R / Python** | Visualization (PCA plots, admixture bars, SNP density) |

---

## Part 1: Principal Component Analysis (PCA)

**What is PCA?**

Principal Component Analysis (PCA) is a statistical method that reduces the dimensionality of complex datasets while retaining most of the variation. In population genetics, PCA is used to visualize genetic relationships between individuals and identify population structure. Individuals that are genetically similar will cluster together in the PCA plot, while genetically distinct populations will separate.

We will use **PLINK** to perform PCA.

### 2.1: Run PCA with PLINK

**Create a file named `run_pca.sh` using `vi`:**

```bash
#!/bin/bash
#SBATCH --job-name=run_pca
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -o run_pca.out
#SBATCH -e run_pca.err

# Load necessary modules
module load plink

# Define input PLINK files base name (pruned data)
INPUT_BASE_PRUNED="my_data_pruned"

echo "Running PCA on ${INPUT_BASE_PRUNED}..."
plink --bfile $INPUT_BASE_PRUNED \
      --pca \
      --out pca_results

echo "PCA complete. Results: pca_results.eigenval, pca_results.eigenvec"
```

**Submit the job:**

```bash
sbatch run_pca.sh
```

This will generate two main output files:

*   `pca_results.eigenval`: Contains the eigenvalues, representing the amount of variance explained by each principal component.
*   `pca_results.eigenvec`: Contains the eigenvectors, which are the principal component scores for each individual.

**❓ Question:** What do the first few principal components (e.g., PC1 and PC2) typically represent in population genetics PCA? How can you use the `eigenval` file to determine how many PCs are important?

### 2.2: Visualize PCA Results (Conceptual)

Visualizing PCA results usually involves plotting the first two or three principal components. This is typically done using statistical software like R or Python with libraries like `ggplot2` or `matplotlib`/`seaborn`.

```python
# Conceptual Python code for PCA visualization
# You would need to transfer pca_results.eigenvec to your local machine
# and use a plotting environment (e.g., Jupyter Notebook, RStudio)

# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns

# Load the eigenvec file (space-separated)
# pca_df = pd.read_csv("pca_results.eigenvec", sep="\s+", header=None)
# Assign column names (e.g., FID, IID, PC1, PC2, ...)
# pca_df.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, pca_df.shape[1] - 1)]

# Example: Plot PC1 vs PC2
# plt.figure(figsize=(8, 6))
# sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="<Population_Group_Column>", s=50) # Replace <Population_Group_Column> if you have population labels
# plt.title("PCA of Genetic Data")
# plt.xlabel(f"PC1 ({pca_df.iloc[0, 2]:.2f}% variance)") # Example for adding variance explained
# plt.ylabel(f"PC2 ({pca_df.iloc[0, 3]:.2f}% variance)")
# plt.grid(True)
# plt.show()
```

**❓ Question:** How would you incorporate known population labels into your PCA plot to validate your findings? What would you look for in the plot to indicate strong population structure?

---

## Part 2: Admixture Analysis

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

