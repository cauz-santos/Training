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

By detecting ROH in our GWAS dataset from Day 7, we can explore how much homozygosity exists in different lines, and which genomic regions may be under breeding pressure.  

### Step 1 — Run bcftools roh

We first detect Runs of Homozygosity (ROH) using the `roh` plugin from **bcftools**.  

```bash
#!/bin/bash
#SBATCH --job-name=roh_bcftools
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH -o roh_bcftools.out
#SBATCH -e roh_bcftools.err

module load bcftools/1.17

VCF="gwas_data_qc.vcf.gz"      # Input from Day 7
OUT="roh_results.txt"

bcftools roh -G30 --rec-rate 1.4e-8 ${VCF} > ${OUT}

echo "ROH calling finished: ${OUT}"
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
grep "RG" roh_results.txt > roh_RG.txt
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

  \[
  F_{ROH} = \frac{\text{SROH}}{\text{Total genome length analyzed}}
  \]

  It represents the **proportion of the genome that is autozygous** (identical by descent).  
  - For example, if an individual has 200 Mb in ROHs and the genome analyzed is 1,000 Mb, then FROH = 0.2 (20%).  
  - Higher FROH indicates more inbreeding.  
  - FROH is comparable across datasets and is often used as a genomic estimate of inbreeding, complementing pedigree-based coefficients.

**Why do we care in breeding?**  
- **High NROH + High FROH** → indicates strong inbreeding, possibly reducing fitness and adaptability.  
- **Low FROH** → indicates more diverse individuals, valuable for crossing and maintaining genetic diversity.  
- These metrics guide breeders in **parental selection, avoiding inbred lines, and managing long-term diversity**.


Open Rstudio and run:
```bash
# ============================
# ROH Analysis in R
# ============================
library(tidyverse)

# Load ROH data (Sample, Chromosome, Length)
roh <- read_delim("roh_RG.txt", delim="\t", skip=1,
                  col_names=c("Sample","Chromosome","Start","End","Length","Other"))

# --- Compute NROH (number of ROHs per individual)
nroh <- roh %>%
  group_by(Sample) %>%
  summarise(NROH = n())

# --- Compute SROH (sum length of ROHs per individual)
sroh <- roh %>%
  group_by(Sample) %>%
  summarise(SROH = sum(Length))

# --- Merge and calculate FROH (fraction of genome in ROH)
genome_size <- 1.2e9   # example genome length; replace with your species
summary_data <- inner_join(nroh, sroh, by="Sample") %>%
  mutate(FROH = SROH / genome_size)

head(summary_data)
```

### Step 4 — Visualizations:  

**A) Scatterplot: NROH vs. SROH**  
In Rstudio:  
```bash
ggplot(summary_data, aes(x=SROH/1e6, y=NROH)) +
  geom_point(color="blue", size=3) +
  theme_minimal() +
  labs(x="Total ROH Length (Mb)", y="Number of ROHs",
       title="Relationship between ROH number and total length")
```

> **Interpretation:**  
> Many short ROHs → older, background inbreeding.  
> Few long ROHs → recent parental relatedness.

**B) Stacked Barplot: ROH Length Categories**  
In Rstudio:  
```bash
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

ggplot(summed_roh, aes(x=Sample, y=total_length/1e6, fill=Category)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(x="Sample", y="Total ROH Length (Mb)", fill="ROH Category",
       title="Distribution of ROH by Length Class") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
```

> **Interpretation:**  
> Short ROHs → reflect background relatedness (drift).  
> Long ROHs → evidence of recent inbreeding.

**C) Boxplot: FROH by Population**    
If you have population metadata (e.g., Sample → Population), you can merge and plot:  

In Rstudio:  
```bash
# Example population file: popmap.txt (Sample \t Population)
popmap <- read_delim("popmap.txt", delim="\t", col_names=c("Sample","Population"))

data_with_pop <- left_join(summary_data, popmap, by="Sample")

ggplot(data_with_pop, aes(x=Population, y=FROH, color=Population)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2, size=2) +
  theme_minimal() +
  labs(y="FROH (Inbreeding Coefficient)", title="Genomic Inbreeding by Population")
```

> **Interpretation:**  
> Higher FROH → stronger inbreeding.  
> Compare populations to identify inbred vs diverse groups.

---

## Session 2 — Linkage Mapping Basics

Here we introduce linkage mapping with both a **toy dataset** and the **Day 7 SNP dataset**.

### Why linkage maps matter

- Provide the **genetic framework** for QTL mapping.  
- Show **marker order and distances** (cM).  
- Detect **recombination patterns** across breeding populations.

### Step 1 — Simulate a toy dataset in R

```r
library(qtl)

# Simulate a small F2 cross with 100 individuals and 3 chromosomes
fake <- sim.cross(map = sim.map(len=c(100,80,60), n.mar=11, anchor=TRUE), n.ind=100, type="f2")

# Estimate recombination fractions
fake <- est.rf(fake)

png("Linkage_toy_rf.png",1000,800)
plotRF(fake)
dev.off()

write.cross(fake, format="csv", filestem="toy_cross")
```

> This generates a small dataset (`toy_cross.csv`) that students can inspect.

### Step 2 — Apply linkage analysis to real SNP data

Convert Day 7 PLINK data to a format suitable for **R/qtl** or **MSTmap**.

```bash
plink --bfile gwas_data_qc --recodeA --out gwas_for_linkage
```

Now load in R:

```r
geno <- read.csv("gwas_for_linkage.raw", sep=" ", header=TRUE)
# Inspect genotypes (0,1,2 dosage coding)
head(geno)
```

> **Note:** Real linkage mapping usually requires a **controlled cross pedigree**, but here we showcase how SNP data can be prepped and visualized for linkage-like patterns.

---

## Session 3 — Machine Learning for Bioinformatics

We’ll train a simple **classification model** on biological data.  
Dataset: We will reuse **phenotypes.csv** (SUC trait) and SNP PCs from Day 6.

### Why ML matters

- ML can integrate **genomic + phenotypic features** to predict traits.  
- Goes beyond linear regression to capture **non-linear patterns**.  
- Used in **genomic prediction** and **trait classification**.

### Step 1 — Prepare dataset

```r
pcs <- read.table("pca_results.eigenvec", header=FALSE)
pheno <- read.csv("phenotypes.csv")
colnames(pcs) <- c("FID","IID",paste0("PC",1:10))
d <- merge(pcs, pheno, by="IID")

head(d)
```

### Step 2 — Train/test split and ML model

We use **random forest** for trait classification (high vs low sucrose).

```r
library(randomForest)

# Define high vs low classes by median split
d$class <- ifelse(d$SUC > median(d$SUC, na.rm=TRUE),"High","Low")

set.seed(42)
train_idx <- sample(1:nrow(d), 0.7*nrow(d))
train <- d[train_idx,]
test <- d[-train_idx,]

rf <- randomForest(as.factor(class) ~ PC1+PC2+PC3+PC4+PC5, data=train, ntree=500)
pred <- predict(rf, test)

table(pred, test$class)
```

### Step 3 — Evaluate accuracy

```r
acc <- mean(pred==test$class)
cat("Test accuracy:", acc, "
")
```

> **Interpretation questions:**  
> - How accurate is the classification?  
> - Which PCs (genomic patterns) seem most predictive?  

> **Relevance:** Even with simple PCs, ML can stratify high vs low sucrose lines. Breeders can expand this with **full SNP data** or **multi-trait models**.

---

# Wrap-up Discussion

- ROH: what does excess homozygosity reveal about breeding history?  
- Linkage: how do recombination patterns constrain breeding?  
- ML: how could predictive models be integrated into selection pipelines?  

---

# Useful Tutorials and Resources

- [PLINK ROH Documentation](https://www.cog-genomics.org/plink/1.9/roh)  
- [qtl R package vignette](https://rqtl.org/tutorials/)  
- [MSTmap for linkage mapping](http://mstmap.org/)  
- [Random Forest in R](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf)  
- [Introduction to Machine Learning for Biologists (YouTube)](https://www.youtube.com/watch?v=tNa99PG8hR8)  
