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

### Step 1 — Call ROH using PLINK

```bash
vi 01_run_roh.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=roh_call
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH -o roh_call.out
#SBATCH -e roh_call.err

module load plink

plink --bfile gwas_data_qc       --homozyg       --homozyg-window-snp 50       --homozyg-snp 50       --homozyg-kb 1000       --homozyg-density 50       --homozyg-gap 100       --out roh_results

echo "ROH calling complete: roh_results.hom"
```

Submit:

```bash
sbatch 01_run_roh.sh
```

- `--homozyg` → tells PLINK to detect runs of homozygosity.  
- `--homozyg-window-snp 50` → minimum number of SNPs in a window.  
- `--homozyg-snp 50` → minimum SNPs required to define a ROH.  
- `--homozyg-kb 1000` → minimum length of ROH (in kb).  
- `--homozyg-density 50` → maximum average kb distance per SNP allowed.  
- `--homozyg-gap 100` → maximum gap between two consecutive homozygous SNPs within a run.  

> Note: The thresholds are adjustable — for plants with high marker density (like GBS or WGS SNPs), you might want stricter settings (e.g., 100 SNPs, 2000 kb). For sparser data (like arrays), you can relax them a bit.

### Step 2 — Summarize ROH

```bash
awk '{print $2,$3,$4,$5}' roh_results.hom | head
```

This file contains: Individual, chromosome, start, end, length, number of SNPs.  

We can summarize by individual:

```bash
awk '{len[$2]+=$6} END{for(i in len) print i,len[i]}' roh_results.hom > roh_per_individual.txt
```

### Step 3 — Visualize in R

```r
roh <- read.table("roh_per_individual.txt")
colnames(roh) <- c("IID","ROH_length")
png("ROH_distribution.png",800,600)
hist(roh$ROH_length/1e6, main="Total ROH length per individual", xlab="ROH (Mb)", breaks=20)
dev.off()
```

> **Interpretation questions:**  
> - Do some individuals carry **longer ROH tracts**?  
> - Could this reflect **inbreeding** or **selection sweeps**?  

> **Relevance:** Breeders use ROH to flag elite lines with excessive autozygosity or to track fixed favorable haplotypes.

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
