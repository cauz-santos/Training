
# Genome-wide Homozygosity Landscape for SNP Selection in Plant Breeding

**Objective:** Visualize genome-wide homozygosity distribution using SNP data to guide marker selection for a DArTseq panel.


## Biological Rationale

In plant breeding, **homozygosity** patterns across the genome can reveal:

- **Fixed regions** (low diversity) → may carry selected alleles or identity-by-descent
- **Variable regions** (higher heterozygosity) → useful for detecting polymorphisms in DArTseq panels
- Avoiding SNPs in **uniformly homozygous** or **paralog-rich** regions improves **panel informativeness**


### Input Data

- Genotype data: `data_pruned.bed/.bim/.fam`
- Reference genome: `EG5_reference_genomic.fna`
- Genome index: `EG5_reference_genomic.fna.fai`
- Working directory: `/lisc/scratch/course/pgbiow/GWAS/plink`


# 🧬 Genome-wide Homozygosity Landscape for SNP Selection in Plant Breeding

**Author:** Luiz A. C. dos Santos  
**Context:** Training in bioinformatics for plant breeding  
**Objective:** Visualize genome-wide homozygosity distribution using SNP data to guide marker selection for a DArTseq panel.

---

### Step 0: Prepare Your Working Directory

Each student should create their own folder for this activity:

```bash
mkdir -p 08_genome_selection/plink
cd 08_genome_selection/plink
```

Then copy all required genotype files from the shared directory:

```bash
cp /lisc/scratch/course/pgbiow/GWAS/plink/data_pruned.* .
```

This ensures they are working with:
- `data_pruned.bed`
- `data_pruned.bim`
- `data_pruned.fam`


### Step 1: Compute Genotype Counts per SNP using PLINK

Create the following SLURM script:

### `01_compute_snp_counts.sh`
```bash
#!/bin/bash
#SBATCH --job-name=snpcounts
#SBATCH --output=snpcounts.out
#SBATCH --error=snpcounts.err
#SBATCH --time=00:05:00
#SBATCH --mem=4G

module load plink

plink \
  --bfile data_pruned \
  --freqx \
  --out snp_counts
```

Submit it:
```bash
sbatch 01_compute_snp_counts.sh
```

This generates the file `snp_counts.frqx` with genotype frequencies across all individuals.


### Step 2: Generate Sliding Windows Across the Genome

We use the `.fai` file of your genome to define 1 Mb windows with 500 kb overlap:

```bash
bedtools makewindows \
  -g /lisc/scratch/course/pgbiow/data/genomes/EG5_reference/EG5_reference_genomic.fna.fai \
  -w 1000000 -s 500000 > genome_1Mb_windows.bed
```

---

## ✅ Step 3: Convert SNP Counts to BED Format

Extract relevant columns from PLINK output:

```bash
awk 'NR>1 {print $1, $4-1, $4, $6, $7, $8}' OFS="\t" snp_counts.frqx > snp_hom.bed
```

Columns:
- $1 = Chromosome
- $4 = SNP position
- $6 = Homozygous REF
- $7 = Heterozygous
- $8 = Homozygous ALT


### Step 4: Aggregate Homozygosity per Window

```bash
bedtools map \
  -a genome_1Mb_windows.bed \
  -b snp_hom.bed \
  -c 4,5,6 \
  -o sum,sum,sum > window_hom_counts.bed
```

This computes the total number of homozygous and heterozygous genotypes in each window.


### Step 5: Visualize Homozygosity in Rstudio

```r
library(readr)
library(dplyr)
library(ggplot2)

# Load BED output
df <- read_tsv("window_hom_counts.bed", col_names = c("chr", "start", "end", "hom_ref", "het", "hom_alt"))

df <- df %>%
  mutate(across(hom_ref:hom_alt, as.numeric)) %>%
  mutate(
    total = hom_ref + het + hom_alt,
    homozygosity = (hom_ref + hom_alt) / total,
    midpoint = (start + end) / 2
  ) %>%
  filter(total > 0)

# Plot per chromosome
ggplot(df, aes(x = midpoint, y = homozygosity)) +
  geom_line(color = "darkblue") +
  facet_wrap(~ chr, scales = "free_x", ncol = 4) +
  labs(title = "Genome-wide Homozygosity Landscape",
       x = "Genomic Position", y = "Proportion Homozygous") +
  theme_minimal()
```


### Interpretation in Plant Breeding

- **Peaks of homozygosity** may indicate:
  - Historical selection (selective sweeps)
  - Autozygosity (inbreeding)
  - Recombination deserts

- **Troughs of homozygosity** often mark:
  - Highly polymorphic regions
  - Good targets for SNP panel design

### Key Takeaways:
- Choose SNPs from regions with **moderate homozygosity** (not fully fixed)
- Avoid markers that fall in **highly homozygous blocks**
- This ensures better resolution and avoids redundant markers in DArTseq panels

---

### Final Outputs

| File                         | Description                              |
|------------------------------|------------------------------------------|
| `snp_counts.frqx`            | Genotype counts for all SNPs             |
| `snp_hom.bed`                | Per-SNP BED-formatted genotype counts    |
| `genome_1Mb_windows.bed`     | Sliding windows across genome            |
| `window_hom_counts.bed`      | Homozygosity counts per window           |
| `02_plot_homozygosity.R`     | R script for genome-wide visualization   |
| `homozygosity_plot.pdf`      | Final plot for teaching and selection    |


### Optional Extensions

- Compare homozygosity across **populations**
- Combine with **LD decay plots** or **GWAS hits**
- Link to **QTL regions** for trait-targeted panel design


