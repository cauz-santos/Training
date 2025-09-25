## Week 2: Applied Bioinformatics for Genomics and Breeding

## Day 8: Runs of Homozygosity (ROH), Linkage Mapping, and Machine Learning

Welcome to **Day 8**!

In this session, we will cover three interconnected topics relevant to plant breeding and applied bioinformatics:

- **Session 1:** Runs of Homozygosity (ROH) from GWAS VCFs (inbreeding, bottlenecks, autozygosity)  
- **Session 2:** Linkage mapping and SNP ordering — toy and real datasets  
- **Session 3:** Machine Learning (ML) basics for biological data, with hands-on practice

The goal is to understand how genomic data can reveal **hidden patterns of inbreeding**, **linkage between loci**, and how **ML models** can be applied to extract predictive power from biological features.

> **Relevance for breeding:**  
> - ROH highlights regions of reduced diversity, important for avoiding inbreeding depression and identifying autozygous blocks carrying key traits.  
> - Linkage mapping is foundational for tracking trait loci and building genetic maps for breeding crosses.  
> - Machine learning provides scalable tools to predict traits, classify samples, and detect complex patterns beyond linear models.  

---

## Session 1 — Runs of Homozygosity (ROH)

Runs of Homozygosity (**ROH**) are long stretches of DNA where both chromosomes are identical (homozygous) across many SNPs. They occur when parents share common ancestors, leading to inbreeding, or when certain regions of the genome are under strong selection (e.g., breeding for specific traits).  

In plant breeding, **ROH are important because**:  
- They reveal **inbreeding levels** within individuals or populations.  
- They help identify **regions under selection**, where breeders may have fixed favorable alleles.  
- They provide insights into **genetic diversity** and **effective population size**.  

By detecting ROH in our GWAS dataset from Day 8, we can explore how much homozygosity exists in different lines, and which genomic regions may be under breeding pressure. 


First create a folder for the file outputs of day 7 in your home directory:
   ```bash
   mkdir 08_rho_linkage
   ```
Now, enter the folder with the command `cd`

Then create some subfolders, for each specific analysis we will perform:
   ```bash
   mkdir vcf roh_runs linkage ml_forest
   ```

### Step 0 — Make sure the VCF exists

**Tip:** ROH works best with dense marker sets. In real analyses you would use the full, unpruned dataset (~600k SNPs).  
👉 For this *practice session*, we will use the file:
`dataset120_chr18.vcf.gz` which is in `/lisc/data/scratch/course/pgbiow/data/VCF`


### Step 1 — Run bcftools roh

We first detect Runs of Homozygosity (ROH) using the `roh` plugin from **bcftools**.  
First, use `cd` to move to directory `08_rho_linkage`

```bash
vi bcftools_roh_job.sh
```

Remember: When the file opens, you are in command mode.
To start typing, press the `i` key (this puts you in insert mode).

Now copy and paste the following script into the file:

```bash
#!/bin/bash
#SBATCH --job-name=roh_bcftools
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH -o roh_bcftools.out
#SBATCH -e roh_bcftools.err

module load bcftools

# Input and output paths
VCF="/lisc/data/scratch/course/pgbiow/data/VCF/dataset120_chr18.vcf.gz"      # Input from Day 7
OUT="roh_runs/roh_results.txt"

# Run ROH calling, estimate AF from all samples
echo "Running bcftools roh..."
bcftools roh -G30 --rec-rate 1.4e-8 --estimate-AF - ${VCF} > ${OUT}

echo "ROH calling finished: ${OUT}"
```

Save and exit (`ESC`, then type `:wq` and press `ENTER`).

Submit the job typing:

```bash
sbatch bcftools_roh_job.sh
```

**Explanation of Parameters:**  
`bcftools roh` → calls the roh plugin to scan the VCF and identify runs of homozygosity.  
`-G30` → sets a genotype quality filter (Phred score ≥ 30). Low-quality genotypes are treated as missing to avoid false ROHs caused by errors.  
`--rec-rate 1.4e-8` → specifies the assumed recombination rate per base per generation. This helps the algorithm decide whether nearby homozygous markers are part of the same ROH or split by recombination.  
*1.4e-8 is a typical genome-wide rate in plants/animals; adjust if species-specific data exists.*  
`${VCF}` → input VCF file.  
`>${OUT}` → saves the full report to roh_results.txt.  

### Step 2 — Extract ROH Segments
The bcftools output contains several record types. ROH blocks are labeled with RG.
We will extract only those:

```bash
# Keep only lines whose first field is exactly "RG" (bcftools ROH segments)
awk -F'\t' '$1=="RG"{print $2,$3,$4,$5,$6}' OFS='\t' ./roh_runs/roh_results.txt > ./roh_runs/roh_segments.tsv
```
This keeps only homozygous blocks with information such as sample ID, chromosome, start, end, length.

### Step 3 — Summarize ROH Statistics in R**  
We’ll now compute and visualize three summary statistics:  

`NROH` → number of ROHs per individual.  
`SROH` → sum of ROH lengths per individual.  
`FROH` → genomic inbreeding coefficient = SROH / genome size.  

**Understanding ROH Statistics**

Once we have detected Runs of Homozygosity (ROHs), we can summarize them with three key statistics that describe the genetic background of each individual:

- **NROH (Number of ROHs per individual)**  
  This is simply how many separate ROH segments are found in a genome.  
  - A **high NROH with mostly short segments** usually reflects background relatedness accumulated over many generations (genetic drift).  
  - A **low NROH but with long segments** suggests recent parental relatedness (e.g., close inbreeding).

- **SROH (Sum of ROH lengths per individual)**  
  This is the total genomic length (in base pairs) covered by ROHs in an individual.  
  - Higher SROH means a larger portion of the genome is homozygous.  
  - Useful to compare overall homozygosity between individuals or populations.

- **FROH (Genomic Inbreeding Coefficient)**  
  This is calculated as:  

 `FROH = SROH / Total genome length analyzed`

  It represents the **proportion of the genome that is autozygous** (identical by descent).  
  - For example, if an individual has 200 Mb in ROHs and the genome analyzed is 1,000 Mb, then FROH = 0.2 (20%).  
  - Higher FROH indicates more inbreeding.  
  - FROH is comparable across datasets and is often used as a genomic estimate of inbreeding, complementing pedigree-based coefficients.

**Why do we care in breeding?**  
- **High NROH + High FROH** → indicates strong inbreeding, possibly reducing fitness and adaptability.  
- **Low FROH** → indicates more diverse individuals, valuable for crossing and maintaining genetic diversity.  
- These metrics guide breeders in **parental selection, avoiding inbred lines, and managing long-term diversity**.

> **To get genome length for FROH**  
> A) Make (or refresh) the FASTA index
>If you don’t already have reference.fa.fai:
> ```bash
> module load samtools
> samtools faidx /lisc/data/scratch/course/pgbiow/data/genomes/date_palm_genome.fna
> ```
> B) Sum the lengths of all contigs in the index
>Total genome length (all contigs)
> ```bash
> awk '{sum+=$2} END{print sum}' /lisc/data/scratch/course/pgbiow/data/genomes/date_palm_genome.fna.fai > genome_length.txt
> cat genome_length.txt
>```

Open Rstudio and run:
```r
# ============================
# ROH Analysis in R
# ============================
library(tidyverse)

# Load ROH data (Sample, Chromosome, Length)
roh <- readr::read_tsv("roh_runs/roh_segments.tsv",
                       col_names=c("Sample","Chromosome","Start","End","Length"),
                       col_types="cciii")

# --- Compute NROH (number of ROHs per individual)
nroh <- roh %>%
  group_by(Sample) %>%
  summarise(NROH = n())

# --- Compute SROH (sum length of ROHs per individual)
sroh <- roh %>%
  group_by(Sample) %>%
  summarise(SROH = sum(Length))

# --- Merge and calculate FROH (fraction of genome in ROH)
genome_size <- 773189301   # example genome length; replace with your species
summary_data <- inner_join(nroh, sroh, by="Sample") %>%
  mutate(FROH = SROH / genome_size)

head(summary_data)
```

### Step 4 — Visualizations:  

**A) Scatterplot: NROH vs. SROH**  
In Rstudio:  
```r
p1 <- ggplot(summary_data, aes(x=SROH/1e6, y=NROH)) +
  geom_point(color="blue", size=3) +
  theme_minimal() +
  labs(x="Total ROH Length (Mb)", y="Number of ROHs",
       title="Relationship between ROH number and total length")

print(p1)  # show in RStudio
```

> **Interpretation:**  
> Many short ROHs → older, background inbreeding.  
> Few long ROHs → recent parental relatedness.

**B) Stacked Barplot: ROH Length Categories**  
In Rstudio:  
```r
# --- ROH categories ---
roh_cat <- roh %>%
  mutate(Category = case_when(
    Length >= 1e6 & Length < 3e6 ~ "1–3 Mb",
    Length >= 3e6 & Length < 5e6 ~ "3–5 Mb",
    Length >= 5e6                ~ ">5 Mb",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Category))

summed_roh <- roh_cat %>%
  group_by(Sample, Category) %>%
  summarise(total_length = sum(Length), .groups="drop")

# --- Plot 2: Distribution by ROH category ---
p2 <- ggplot(summed_roh, aes(x=Sample, y=total_length/1e6, fill=Category)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(x="Sample", y="Total ROH Length (Mb)", fill="ROH Category",
       title="Distribution of ROH by Length Class") +
  theme(axis.text.x=element_text(angle=45,hjust=1))

print(p2)
```

> **Interpretation:**  
> Short ROHs → reflect background relatedness (drift).  
> Long ROHs → evidence of recent inbreeding.
>
> Short/medium ROHs reflect the base germplasm diversity.
> Long ROHs are often introduced intentionally during line development.

**C) Boxplot: FROH by Population**    
If you have population metadata (e.g., Sample → Population), you can merge and plot:  

In Rstudio:  
```r
# --- Load population table and merge ---
popmap <- read_csv("/lisc/data/scratch/course/pgbiow/data/metadata/gwas_pop_table_120.csv")

data_with_pop <- summary_data %>%
  left_join(popmap, by = c("Sample" = "IID"))

# --- Plot 3: FROH per population ---
p3 <- ggplot(data_with_pop, aes(x = Population, y = FROH, color = Population)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  theme_minimal(base_size = 14) +
  labs(y = "FROH (Inbreeding Coefficient)", 
       title = "Genomic Inbreeding by Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)

# --- Save all plots as PDF in roh_runs/ ---
pdf("roh_runs/roh_plots.pdf", width=10, height=7)
print(p1)
print(p2)
print(p3)
dev.off()

cat("All plots saved to roh_runs/roh_plots.pdf\n")
```

> **Interpretation:**  
> Higher FROH → stronger inbreeding.  
> Compare populations to identify inbred vs diverse groups.
>
> Populations with higher FROH have reduced genetic diversity, which can lead to inbreeding depression (loss of vigor, yield).

**D) Outliers**  
Now we will go beyond population averages and identify **specific individuals** that deviate strongly from the expected inbreeding levels.  

- Individuals with **very high FROH** are likely the result of **close-relative mating or selfing**, and may carry risks of **inbreeding depression**.  
- Individuals with **very low FROH** could indicate **sample contamination, mislabeling, or unexpected admixture**.  
- We will also flag individuals with an **excess of very long ROHs**, which are signatures of **recent inbreeding events**.  

```r
# ============================
# Identify Outlier Individuals
# ============================

# --- Step 1: Calculate thresholds for FROH ---
froh_mean <- mean(summary_data$FROH, na.rm = TRUE)
froh_sd   <- sd(summary_data$FROH, na.rm = TRUE)

high_cutoff <- froh_mean + 2 * froh_sd
low_cutoff  <- froh_mean - 2 * froh_sd

cat("FROH thresholds:\n",
    "High outlier >", round(high_cutoff, 3), "\n",
    "Low outlier <", round(low_cutoff, 3), "\n")

# --- Step 2: Flag outliers ---
outliers <- summary_data %>%
  mutate(
    OutlierReason = case_when(
      FROH > high_cutoff ~ "High FROH (possible inbreeding)",
      FROH < low_cutoff ~ "Low FROH (possible contamination/mislabel)",
      SROH > quantile(SROH, 0.95, na.rm = TRUE) ~ "Excess total ROH length",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(OutlierReason)) %>%
  left_join(popmap, by = c("Sample" = "IID"))

# --- Step 3: Save outlier table ---
readr::write_csv(outliers, "roh_runs/outliers.csv")

# --- Step 4: Print preview for students ---
print(outliers)

cat("\n Outlier individuals saved to roh_runs/outliers.csv\n")
```

**Discussion Questions:**  
- What does a high NROH but low SROH suggest about a population’s history?
- Why do very long ROHs point to recent parental relatedness?
- How does FROH help breeders avoid inbreeding depression?
- How could breeders use ROH patterns to balance fixing favorable traits vs maintaining diversity?

**Relevance to Breeding:**  
- **Decision-making:** Avoid crossing individuals with high FROH to reduce inbreeding depression.  
- **Marker discovery:** Long ROHs overlapping trait regions may indicate selected alleles.
- **Conservation:** Populations with lower ROH diversity are more resilient.
- **Breeding strategy:** Helps decide when to outcross vs self/cross inbred lines.

---

## Session 2 — Linkage Mapping Basics
A **linkage map** is a genetic map that shows the relative positions of markers (e.g., SNPs, microsatellites) along chromosomes, based on how often they are inherited together.  
The unit is the **centimorgan (cM)**, which represents the recombination frequency:  
- **1 cM ≈ 1% chance** that a crossover occurs between two markers during meiosis.  

Because recombination is less likely between markers that are close together, linkage maps help us order markers and estimate distances between them.  

### Why linkage maps are important in breeding

- Provide the **genetic framework** for locating QTLs (quantitative trait loci).  
- Help breeders track recombination events and **design crosses more effectively**.  
- Allow identification of **linkage drag** (undesirable alleles inherited together with target traits).  
- Serve as the foundation for integrating **genomics with classical breeding**.  

In this session, we’ll practice linkage mapping first with a **toy dataset** (F2 cross simulation) and then explore how to prepare our real SNP dataset (from Day 7) for linkage-like analyses.

### Linkage Map vs LD Map

Although the terms sound similar, **linkage maps** and **LD maps** describe different concepts.

| Feature              | Linkage Map                              | LD Map                              |
|----------------------|-------------------------------------------|--------------------------------------|
| **Population type**  | Controlled cross (F2, RIL, DH, BC)        | Natural/diversity panel              |
| **Measure**          | Recombination frequency (cM)              | Statistical correlation (r², D’)     |
| **Output**           | Marker order + distances (map)            | LD blocks/haplotype structure        |
| **Reflects**         | True meiotic recombination                | Recombination + demography + drift   |
| **Applications**     | QTL mapping, breeding populations         | GWAS, diversity, haplotype tagging   |


### Step 1 — Simulate a toy dataset (F2 population)

Because linkage maps require **pedigreed crosses** (e.g., F2), we’ll simulate one in R to show how maps are built.

```r
## ==========================================
## Linkage Mapping Basics (simulation)
## ==========================================
library(ggplot2)
library(reshape2)

set.seed(123)

## -----------------------
## Parameters
## -----------------------
N_ind           <- 300   # "sample size" controlling noise
n_chr           <- 5     # number of chromosomes
markers_per_chr <- 20    # markers per chromosome
min_gap_cM      <- 2
max_gap_cM      <- 10

## -----------------------
## Haldane functions
## -----------------------
haldane_r_from_cM <- function(d_cM) {
  d_M <- d_cM / 100
  0.5 * (1 - exp(-2 * d_M))
}
haldane_cM_from_r <- function(r) {
  r <- pmin(pmax(r, 0), 0.499999)
  d_M <- -0.5 * log(1 - 2 * r)
  100 * d_M
}

## -----------------------
## 1) Build a true linkage map
## -----------------------
chroms <- paste0("Chr", seq_len(n_chr))
make_chr_map <- function(chr, m, min_gap, max_gap) {
  gaps <- runif(m, min_gap, max_gap)
  cM   <- sort(cumsum(gaps))
  data.frame(Chromosome = chr,
             Marker = paste0(chr, "_M", seq_len(m)),
             cM = cM,
             stringsAsFactors = FALSE)
}
true_map <- do.call(rbind, lapply(chroms, make_chr_map,
                                  m = markers_per_chr,
                                  min_gap = min_gap_cM,
                                  max_gap = max_gap_cM))

markers <- true_map$Marker
n_mark  <- length(markers)
chr_of  <- setNames(true_map$Chromosome, true_map$Marker)
cM_of   <- setNames(true_map$cM, true_map$Marker)

## -----------------------
## 2) Simulate pairwise r with noise
## -----------------------
r_exp <- matrix(0.5, n_mark, n_mark, dimnames = list(markers, markers))
for (i in seq_len(n_mark)) {
  for (j in seq_len(n_mark)) {
    if (i == j) {
      r_exp[i, j] <- 0
    } else if (chr_of[markers[i]] == chr_of[markers[j]]) {
      d_cM <- abs(cM_of[markers[i]] - cM_of[markers[j]])
      r_exp[i, j] <- haldane_r_from_cM(d_cM)
    }
  }
}
noise_sd <- sqrt(r_exp * (1 - r_exp) / N_ind)
r_obs    <- pmin(pmax(r_exp + rnorm(length(r_exp), 0, as.vector(noise_sd)), 0), 0.5)
dimnames(r_obs) <- list(markers, markers)

## -----------------------
## 3) Plot recombination fraction heatmap
## -----------------------
r_melt <- reshape2::melt(r_obs, varnames = c("Marker1", "Marker2"), value.name = "r")

ggplot(r_melt, aes(x = Marker1, y = Marker2, fill = r)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "yellow", name = "Recombination fraction (r)",
                      limits = c(0, 0.5)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5),
        axis.text.y = element_text(size = 5)) +
  labs(title = "Recombination Fraction Heatmap",
       x = "Markers", y = "Markers")

dir.create("linkage", showWarnings = FALSE)
ggsave("linkage/r_heatmap.pdf", width = 9, height = 9)

# ============================
# Linkage Map Plot 
# ============================

# Range for y-axis (genetic distances in cM)
yl <- range(0, tapply(true_map$cM, true_map$Chromosome, max))

# Preview in RStudio
plot(NA, xlim = c(0.5, n_chr + 0.5), ylim = yl,
     xlab = "Chromosome", ylab = "Genetic distance (cM)",
     xaxt = "n", main = "Linkage Map (true cM)")
axis(1, at = seq_len(n_chr), labels = chroms)

for (i in seq_along(chroms)) {
  chr <- chroms[i]
  chr_map <- true_map[true_map$Chromosome == chr, ]
  # chromosome vertical line
  lines(c(i, i), c(0, max(chr_map$cM)), lwd = 3, col = "grey40")
  # marker positions
  points(rep(i, nrow(chr_map)), chr_map$cM, pch = 19, col = "red")
  # marker labels
  text(rep(i + 0.15, nrow(chr_map)), chr_map$cM,
       labels = chr_map$Marker, cex = 0.7, pos = 4)
}

# ============================
# Save to PDF
# ============================
pdf("linkage/linkage_map_true.pdf", width = 10, height = 7)
plot(NA, xlim = c(0.5, n_chr + 0.5), ylim = yl,
     xlab = "Chromosome", ylab = "Genetic distance (cM)",
     xaxt = "n", main = "Linkage Map (true cM)")
axis(1, at = seq_len(n_chr), labels = chroms)

for (i in seq_along(chroms)) {
  chr <- chroms[i]
  chr_map <- true_map[true_map$Chromosome == chr, ]
  lines(c(i, i), c(0, max(chr_map$cM)), lwd = 3, col = "grey40")
  points(rep(i, nrow(chr_map)), chr_map$cM, pch = 19, col = "red")
  text(rep(i + 0.15, nrow(chr_map)), chr_map$cM,
       labels = chr_map$Marker, cex = 0.7, pos = 4)
}
dev.off()


cat(" Outputs saved in linkage/: r_heatmap.pdf and linkage_map_true.pdf\n")
```

**Outputs:**  
- `r_heatmap.pdf` → heatmap of recombination fractions.  
- `linkage_map_true.pdf` → visualization of chromosomes with marker order and genetic distances (cM), i.e. the constructed linkage map.

> **How to Interpret the RF Heatmap**  
> - **Dark blocks (low rf)** → tightly linked markers.   
> - **Light blocks (high rf)** → unlinked markers (different chromosomes or very distant on the same chromosome).  
> - You’ll see **5 dark diagonal blocks** (one per chromosome), which shows that markers group correctly into **linkage groups**.  
>
> 👉 This mimics what happens in **real linkage mapping**: you cluster markers into groups (chromosomes) based on recombination.


**Relevance for breeding:**  
- Shows how **recombination frequency translates into distances**.  
- Demonstrates that **markers closer together recombine less frequently**.  
- Provides the foundation for **QTL mapping**.


## Step 2 — Apply linkage-like analysis to real SNP data
Our GWAS dataset from Day 7 (`gwas_data_qc.*`) is **not from an F2 cross**, so we cannot build a traditional linkage map.  
But we can still explore **linkage disequilibrium (LD) patterns** in this diversity panel.  
LD reflects **correlations between SNPs** due to shared ancestry, drift, or selection.


**1) Convert PLINK data to allele dosage format**  
We need SNP genotypes coded as **0, 1, 2** (number of alternate alleles).  

```bash
cd linkage
```

Run:
```bash
module load PLINK
```

```bash
plink --bfile /lisc/scratch/course/pgbiow/backup/07_gwas_selection/plink/gwas_data_qc --recodeA --allow-extra-chr --out gwas_for_linkage
```

**Output:**  
- `gwas_for_linkage.raw` → matrix of individuals × SNPs with dosage coding.  
- Columns: `FID`, `IID`, `PAT`, `MAT`, `SEX`, `PHENOTYPE`, followed by SNPs.



**2) Load genotype matrix into R**  
In Rstudio:
```r
geno <- read.csv("./linkage/gwas_for_linkage.raw", sep=" ", header=TRUE)

# Inspect first columns
head(geno[,1:10])
```

**Let’s break it down:**    
`FID / IID` → Family ID and Individual ID (your sample names).  
`PAT / MAT / SEX / PHENOTYPE` → placeholders for pedigree, sex, and phenotype (here all 0 or -9 because you didn’t provide this info).  

The next columns are the SNP markers. Each column corresponds to one SNP (named by chromosome/position/alleles).  
Example: `chrLG1.83953.G.T_T` → SNP on chromosome LG1 at position 83,953 with alleles G/T.  
The values (0, 1, 2, NA) are genotypes coded as allele counts:  
               `0` = homozygous for the reference allele  
               `1` = heterozygous  
               `2` = homozygous for the alternate allele  
               `NA` = missing genotype  

**3) Compute pairwise LD with PLINK**  
PLINK can directly compute LD between SNP pairs:

```bash
plink --bfile /lisc/scratch/course/pgbiow/backup/07_gwas_selection/plink/gwas_data_qc --r2 --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0 --allow-extra-chr --out gwas_ld
```

**Outputs:**  
- `gwas_ld.ld` → file with SNP pairs, distance, and r².


**4) Visualize LD decay in R**  
We can check how LD decays with physical distance:

In Rstudio:
```r
library(data.table)
library(ggplot2)
library(dplyr)

# Load LD file
ld <- fread("./linkage/gwas_ld.ld")

# Subsample for speed
ld_sub <- ld[sample(.N, min(20000, .N))]

# Scatterplot with GAM smoother
p1 <- ggplot(ld_sub, aes(x = BP_B - BP_A, y = R2)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red") +
  theme_minimal(base_size = 14) +
  labs(x = "Distance between SNPs (bp)", y = expression(r^2),
       title = "LD Decay (subsampled SNP pairs)")

# Binned averages (10 kb bins)
ld_binned <- ld %>%
  mutate(bin = cut(BP_B - BP_A, breaks = seq(0, max(BP_B - BP_A), by = 10000))) %>%
  group_by(bin) %>%
  summarise(mean_dist = mean(BP_B - BP_A), mean_r2 = mean(R2, na.rm=TRUE))

p2 <- ggplot(ld_binned, aes(x = mean_dist, y = mean_r2)) +
  geom_point(color = "blue") +
  geom_line(color = "red") +
  theme_minimal(base_size = 14) +
  labs(x = "Distance between SNPs (bp)", y = expression(r^2),
       title = "LD Decay (10 kb binned averages)")

# Save both plots to PDF
pdf("./linkage/LD_decay_plots.pdf", width = 10, height = 6)
print(p1)
print(p2)
dev.off()
```

**Output:**  
- A scatterplot showing how r² decreases as SNP distance increases.  
- The smoother line shows the **LD decay curve**.

> NOTE: The LD decay plot shows how linkage disequilibrium (r²) declines with increasing physical distance between SNPs. In our date palm populations, r² values remain low (≤0.08) even at short distances, indicating very rapid LD decay. This is a hallmark of cross-pollinated crops, where recombination quickly breaks associations between markers. For breeders, this means that mapping resolution is potentially very high, causal genes can be pinpointed with dense SNP data, but it also requires larger marker panels to capture meaningful associations, unlike in self-pollinating crops (e.g., rice or wheat) where LD extends across much longer genomic regions.


**5) LD heatmap for one chromosome region**  
We can zoom in and visualize LD as a heatmap:

In Rstudio:
```r
library(data.table)
library(reshape2)
library(ggplot2)
library(stringr)

# ---- 1) Read a manageable SNP subset (adjust 7:206 as you like) ----
geno <- fread("./linkage/gwas_for_linkage.raw")
snp_subset <- geno[, 7:206, with = FALSE]           # ~200 SNPs is readable
snps <- colnames(snp_subset)

# ---- 2) LD matrix (r^2) ----
ld_mat <- cor(snp_subset, use = "pairwise.complete.obs")^2

# ---- 3) Long format + use numeric indices on axes ----
ld_df <- reshape2::melt(ld_mat, varnames = c("SNP_A","SNP_B"), value.name = "R2")
ld_df$A_i <- as.integer(factor(ld_df$SNP_A, levels = snps))
ld_df$B_i <- as.integer(factor(ld_df$SNP_B, levels = snps))

# Keep upper triangle only (classic LD view). Comment this line if you want full square.
ld_df <- ld_df[ld_df$A_i < ld_df$B_i, ]

# ---- 4) Sparse, truncated tick labels every 10 SNPs ----
idx_breaks <- seq(1, length(snps), by = 10)
idx_labels <- str_trunc(snps[idx_breaks], 18)       # shorten long IDs

p <- ggplot(ld_df, aes(A_i, B_i, fill = R2)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_gradient(limits = c(0,1), low = "white", high = "red", name = expression(r^2)) +
  scale_x_continuous(breaks = idx_breaks, labels = idx_labels, expand = c(0,0)) +
  scale_y_continuous(breaks = idx_breaks, labels = idx_labels, expand = c(0,0)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "LD Heatmap (subset)", x = "SNP (index)", y = "SNP (index)")

print(p)

# Save a clean PDF
ggsave("./linkage/LD_heatmap_subset.pdf", p, width = 6.5, height = 5.5)
```

**Output:**  
- Heatmap where **dark red = strong LD** and **pink/white = weak LD**.  

> **Interpretation:**  
> **LD decay plot:**  
> - Fast decay → high recombination, more diversity.  
> - Slow decay → strong structure, less recombination, possible selection.  
>
> **LD heatmap:**
> - Dark red (close to 1.0) → SNP pairs are in strong LD, meaning they are almost always inherited together (little to no recombination between them in your population).
> - Light pink/white (close to 0.0) → SNP pairs are in weak or no LD, meaning they assort independently (lots of recombination, or just no correlation).

**Key point:**   
Even though we cannot build a **linkage map** without a controlled F2 cross,  
we can still use LD to study **marker correlations**, which is the foundation of **GWAS** and **haplotype-based breeding**.

**Discussion Questions:**

- Why do we need **F2 or RIL populations** to build linkage maps?  
- What information can LD maps provide that linkage maps cannot (and vice versa)?  
- How could you use both maps in a breeding program?  

---

## Session 3 — Basics of Machine Learning for Bioinformatics 
In this session we will explore how **machine learning (ML)** can be applied in bioinformatics and plant breeding.  
We will train a simple **classification model** using genomic principal components (PCs, from Day 6) and sucrose phenotype data (from Day 7).


## Why ML matters in breeding and bioinformatics
- ML can integrate **genomic + phenotypic features** to predict traits.  
- Goes beyond linear regression to capture **non-linear patterns**.  
- Used in **genomic prediction**, **trait classification**, and **selection of top candidates**.


## What is Random Forest?
**Random Forest** is an ensemble machine learning algorithm. It builds **many decision trees** on random subsets of the data and variables. Each tree makes a prediction, and the forest combines them by majority vote (classification) or averaging (regression).  

Why it works well:
- Handles **non-linear relationships** between SNP/genomic data and traits.  
- Reduces **overfitting** because no single tree dominates.  
- Works with **high-dimensional data** (many predictors, like SNP PCs).  

> 👉 Think of it as a “committee of experts”: each decision tree is one expert, and the forest is the combined decision. This makes the predictions more robust than relying on just one.


**Step 0 - Installing the `randomForest` package in RStudio**  

Before running the exercise, make sure the `randomForest` package is installed in your **personal R library**.  
This is necessary because the global R installation on the cluster is read-only, and each user must install their own copy of additional packages.

In RStudio, run the following once:

```r
# Create a personal library path (only needed once)
dir.create("/lisc/home/user/<your_username>/R/x86_64-pc-linux-gnu-library/4.5",
           recursive = TRUE, showWarnings = FALSE)

# Tell R to use your personal library
.libPaths("/lisc/home/user/<your_username>/R/x86_64-pc-linux-gnu-library/4.5")

# Install randomForest into your library
install.packages("randomForest", repos = "https://cloud.r-project.org",
                 lib = "/lisc/home/user/<your_username>/R/x86_64-pc-linux-gnu-library/4.5")

# Load the package to test installation
library(randomForest)
```

**Step 1 — Prepare dataset**  

We’ll merge PCA results (Day 6) with sucrose phenotypes (Day 7).

In Rstudio:
```r
library(randomForest)

# Load data
pcs <- read.table("/lisc/scratch/course/pgbiow/06_diversity_structure/pca/pca_results.eigenvec",
                  header=FALSE)
pheno <- read.csv("/lisc/scratch/course/pgbiow/metadata/gwas_phen_table_120.csv")

# Name PCA columns (first two are IDs, then PCs)
colnames(pcs) <- c("FID","IID",paste0("PC",1:10))

# Merge on IID
d <- merge(pcs, pheno, by="IID")
```

**Step 2 — Clean and define classes:**  
We first remove samples that have missing phenotypes, then classify lines into High vs Low sucrose using a median split.

```r
# Remove missing phenotypes
d_clean <- d[!is.na(d$SUC), ]
cat("Samples after removing missing SUC:", nrow(d_clean), "\n")
cat("Dropped samples:", nrow(d) - nrow(d_clean), "\n")

# Define sucrose class
d_clean$class <- ifelse(d_clean$SUC > median(d_clean$SUC, na.rm=TRUE), "High", "Low")
```

**Step 3 — Train/test split:**  
We split the dataset into 70% training and 30% test samples.

```r
set.seed(42)
train_idx <- sample(1:nrow(d_clean), 0.7*nrow(d_clean))
train <- d_clean[train_idx, ]
test  <- d_clean[-train_idx, ]

cat("Training samples:", nrow(train), "\n")
cat("Testing samples:", nrow(test), "\n")
```

**Step 4 —Train Random Forest:**  
We train a Random Forest classifier using the 10 PCs as predictors.

In Rstudio:
```r
library(randomForest)

rf <- randomForest(as.factor(class) ~ PC1+PC2+PC3+PC4+PC5,
                   data=train, ntree=500, importance=TRUE)

# Predictions on test set
pred <- predict(rf, test)
```

**Step 3 — Evaluate performance:**  
After training a machine learning model, we need to test how well it performs on data it has **never seen before** (the test set).  
One simple way is to calculate the **accuracy**: the proportion of correct predictions compared to the true classes.

In Rstudio:
```r
# Accuracy
acc <- mean(pred == test$class)
cat("Test accuracy:", round(acc, 3), "\n")

# Confusion matrix
cat("Confusion matrix:\n")
print(table(Predicted=pred, Actual=test$class))
```

**What this does:**  
- Compares the model’s predicted sucrose classes (`pred`) to the true classes (`test$class`).  
- Calculates the proportion of correct predictions with `mean(pred == test$class)`.  
- Displays a **confusion matrix**, which shows how many High and Low samples were classified correctly vs. incorrectly.  

**Why this matters:**  
- Accuracy tells us how well the model can generalize to unseen samples.  
- A confusion matrix gives deeper insight: e.g. the model may predict "Low" sucrose better than "High".  
- If accuracy is high (>0.7–0.8), the PCs are informative for distinguishing sucrose classes.  
- If accuracy is low (<0.6), sucrose differences may not be strongly captured by population structure (PCs), or the model may need more features.  

**How to interpret:**  
- Accuracy is a quick, intuitive measure of performance, but it’s not the whole story.  
- A balanced dataset (equal High/Low classes) makes accuracy meaningful; otherwise, accuracy alone can be misleading.  
- In real genomic prediction, we’d also use **cross-validation** for stability and additional metrics:  
  - **AUC** (Area Under Curve) → for classification tasks.  
  - **RMSE** (Root Mean Square Error) → for regression on continuous traits.  
- For training purposes, accuracy + confusion matrix provide a straightforward way to evaluate the model.  


**Step 4 (Optional) — Variable importance and plots**  
In many machine learning models, especially **Random Forests**, we can check **which features (variables)** are most useful for making predictions. This is called **variable importance**.

In Rstudio:
```r
# Variable importance (which PCs are most predictive?)
importance(rf)
varImpPlot(rf)
```

**--- Save evaluation artifacts (Optional) ---**  
It is good practice to save your model evaluation results, so they can be reused or inspected later.
```r
# Confusion matrix + accuracy
cm <- table(Pred = pred, True = test$class)
print(cm)
write.table(as.data.frame(cm), "RF_confusion_matrix.tsv", sep = "\t", row.names = FALSE)
writeLines(sprintf("Test_accuracy\t%.3f", acc), "RF_accuracy.txt")

# Variable importance (table + plot)
vip <- importance(rf)
write.csv(vip, "RF_variable_importance.csv", row.names = TRUE)
png("RF_varImp.png", width = 1000, height = 800, res = 150)
varImpPlot(rf, main = "Random Forest Variable Importance (PCs)")
dev.off()

# Reproducibility breadcrumbs
set.seed(42)
sink("R_sessionInfo.txt"); print(sessionInfo()); sink()
```

**For what purpose:**  
- To identify the most informative genomic patterns (e.g., PC1 vs PC3).  
- To guide breeders or researchers toward which parts of the genome contribute most to trait differences.
- To avoid overfitting by focusing on the strongest predictors.  

**How to interpret:**  
The plot ranks PCs by their contribution to prediction accuracy.  
- If PC1 is most important → the strongest population structure is predictive of sucrose levels.
- If PC3/PC4 are more important → subtle genomic patterns explain sucrose variation.

### Interpretation questions
- How accurate is the classification?  
- Which PCs (genomic patterns) are most predictive?  

### Relevance to Breeding
Even with simple PCs, ML can stratify **high vs low sucrose lines**.  
In practice, breeders can:  
- Use **Random Forest or other ML models** with **full SNP data**.  
- Apply ML for **genomic prediction**, ranking candidates before phenotyping.  
- Integrate **multi-trait and environmental data** for more robust selection.  





---
## Conclusion

Today, you **quantified inbreeding with ROH (NROH/SROH/FROH)**, contrasted **recombination (linkage) with correlation (LD)** to understand block structure, and learned how these signals guide marker spacing and cross design.
You also trained a **lightweight Random Forest on PCs to rank high/low sucrose lines** — an immediately useful triage that can later be upgraded to full genomic prediction.


## You have completed **Day 8**!

---

### Useful Tutorials and Resources

- [Heterozygosity and runs of homozygosity](https://anopheles-genomic-surveillance.github.io/workshop-5/module-4-roh.html)  
- [Detecting runs of homozygosity (RoH)](https://samtools.github.io/bcftools/howtos/roh-calling.html)
- [MSTmap for linkage mapping](http://mstmap.org/)
- [OneMap](https://github.com/augusto-garcia/onemap)  
- [Random Forest in R](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf)  
- [Machine Learning for Biologists - Guide](https://www.nature.com/articles/s41580-021-00407-0)
- [Machine Learning for Biologists](https://carpentries-incubator.github.io/ml4bio-workshop/01-introduction/index.html)



