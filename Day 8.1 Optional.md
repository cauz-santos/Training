# Genomic Selection and Heterosis Prediction

## Introduction

In this practical, we will perform **genomic selection (GS)** and **heterosis prediction** using the provided Verdant datasets of parents and offspring.  
The goal is to learn how to train GS models on **offspring SNP + phenotype data**, and then use these models to predict the **best parent combinations** that maximize heterosis.

Each participant should run this practical inside their own folder:

Please create the folder for:

```
mkdir 08_genomic_selection/
```


### Step 1 - Prepare the Environment

Load required R packages:

```{r setup, message=FALSE}
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(rrBLUP)
library(purrr)
```

Define the working directory (adapt if needed):

```{r paths}
base <- "/lisc/scratch/course/pgbiow/data/GS_Heterosis/Heterosis_prediction_verdant"
out_dir <- file.path(base, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
```


### Step 2 - Analysis

We will use:

- **Parent genotypes** from: `training_sets/parents/parents.raw`
- **Offspring genotypes** from: `training_sets/offsprings/offsprings.raw`
- **Offspring phenotypes** from: `testing_sets/list_input_parent_pairs_offsprings_*.csv`
- **Parent phenotypes** from: `training_sets/ind_names_phenotypes_heterosis_Sf_30.csv`

```r
#!/usr/bin/env Rscript

# ==============================
# Day 08: Genomic Selection & Heterosis Prediction
# ==============================

# Load libraries from user space
.libPaths(Sys.getenv("R_LIBS_USER"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rrBLUP)
  library(caret)
})

# ------------------------------
# Paths (adapt to cluster)
# ------------------------------
base <- "/lisc/scratch/course/pgbiow/data/GS_Heterosis/Heterosis_prediction_verdant"

parent_geno  <- file.path(base, "training_sets/parents/Report_DOp24-9710_SNP.csv")
offspring_geno <- file.path(base, "training_sets/offsprings/Report_DOp25-10424_SNP.csv")
offspring_pheno <- file.path(base, "testing_sets/list_input_parent_pairs_offsprings_DOp25-10620.csv")
parent_pheno <- file.path(base, "training_sets/ind_names_phenotypes_heterosis_Sf_30.csv")

out_dir <- file.path(base, "results_day08")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TRAIT <- "MBW-MA"   # change to target other traits

# ------------------------------
# Load data
# ------------------------------
cat("Reading genotype data...\n")
parents <- fread(parent_geno)
offspring <- fread(offspring_geno)

cat("Reading phenotypes...\n")
offs_ph <- read_csv(offspring_pheno, show_col_types = FALSE)
pars_ph <- read_csv(parent_pheno, show_col_types = FALSE)

# Match SNPs between parents and offspring
snps <- intersect(parents$ID, offspring$ID)
parents <- parents[ID %in% snps]
offspring <- offspring[ID %in% snps]

cat("Common SNPs: ", length(snps), "\n")

# ------------------------------
# Build offspring training set
# ------------------------------
offs_geno <- as.matrix(offspring[, -1, with = FALSE])
rownames(offs_geno) <- offspring$ID

# Select trait
offs_pheno <- offs_ph %>%
  transmute(ID = as.character(id),
            Trait = as.numeric(.data[[TRAIT]])) %>%
  filter(!is.na(Trait))

# Align
idx <- match(offs_pheno$ID, rownames(offs_geno))
offs_geno <- offs_geno[idx, , drop=FALSE]

y <- offs_pheno$Trait
Z <- scale(offs_geno, center = TRUE, scale = TRUE)

# ------------------------------
# Train GS model with rrBLUP
# ------------------------------
cat("Training rrBLUP model...\n")
fit <- mixed.solve(y = y, Z = Z)
beta <- as.vector(fit$u)

# ------------------------------
# Predict crosses
# ------------------------------
cat("Predicting crosses...\n")
par_geno <- as.matrix(parents[, -1, with = FALSE])
rownames(par_geno) <- parents$ID
par_ids <- rownames(par_geno)

combos <- t(combn(par_ids, 2))
results <- data.frame(Parent1 = combos[,1], Parent2 = combos[,2])

preds <- apply(combos, 1, function(x) {
  g1 <- par_geno[x[1],]
  g2 <- par_geno[x[2],]
  gmean <- (g1 + g2) / 2
  sum(gmean * beta, na.rm = TRUE)
})
results$Hybrid_pred <- preds

# ------------------------------
# Compute heterosis
# ------------------------------
id_col <- intersect(names(pars_ph), c("ID","id","name","genotype"))[1]
trait_col <- setdiff(names(pars_ph), c(id_col,"pop","group"))

pars_ph2 <- pars_ph %>%
  transmute(ID = as.character(.data[[id_col]]),
            Trait = as.numeric(.data[[trait_col]]))

id2p <- setNames(pars_ph2$Trait, pars_ph2$ID)

results <- results %>%
  mutate(P1 = id2p[Parent1],
         P2 = id2p[Parent2]) %>%
  filter(!is.na(P1), !is.na(P2)) %>%
  mutate(MP = (P1 + P2)/2,
         BP = pmax(P1, P2),
         MPH = (Hybrid_pred - MP)/MP,
         BPH = (Hybrid_pred - BP)/BP) %>%
  arrange(desc(MPH))

# ------------------------------
# Save results
# ------------------------------
write_csv(results, file.path(out_dir, paste0("GS_Heterosis_", TRAIT, ".csv")))

# Heatmap
library(ggplot2)
p <- ggplot(results[1:100,], aes(x=Parent1, y=Parent2, fill=MPH)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_minimal(base_size=8) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  labs(title=paste0("Predicted MPH (top 100 crosses) - ", TRAIT))

ggsave(file.path(out_dir, paste0("Heatmap_MPH_", TRAIT, ".png")), p, width=12, height=10)

cat("Done! Results saved in ", out_dir, "\n")
```


## Conclusion

You have now:

- Trained a **GS model with additive + dominance effects** on the offspring data  
- Predicted **all possible parent × parent crosses**  
- Computed **heterosis measures (MPH, BPH)**  
- Produced a **ranked CSV** and a **heatmap** for visualization  
