---
# title: "Genomic Selection Practical: Predicting Heterosis & Choosing Parent Combinations (Simulated Data)"

---

## Overview

**Goal:** hands-on practical to (1) train a genomic prediction model (GBLUP) on simulated data, (2) obtain genomic estimated breeding values (GEBVs) for inbred parents, (3) predict hybrid performance for all possible crosses, (4) compute heterosis (MPH and BPH), and (5) select the best parent combinations.

This notebook is **self-contained**: it **simulates** genotypes and phenotypes for 30 inbred parents and walks through the pipeline.

**You’ll produce:**
- A table of **predicted hybrid values** and **heterosis** (MPH and BPH).
- A **heatmap** of predicted heterosis for the full crossing matrix.
- A ranked table of the **top crosses**.

> Default model is **additive GBLUP** via `rrBLUP`. An optional **additive + dominance** section using `sommer` is provided at the end if you want to go one step further.

# Setup

```r
# Install (run once if needed)
# install.packages(c("rrBLUP", "ggplot2", "dplyr", "tidyr", "Matrix", "reshape2"))

library(rrBLUP)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(reshape2)
set.seed(123)
```


# 1) Simulate genotypes for 30 inbred parents

We simulate **M=10000 bi-allelic SNPs** for **N=30** fully-inbred parents (genotypes coded as 0/2 for homozygotes; 1 rarely appears due to inbreeding). Minor allele frequencies (MAF) are sampled to get realistic variation.

```r
N <- 30    # number of parents
M <- 10000  # number of SNPs

# Draw minor allele frequencies
maf <- runif(M, min = 0.05, max = 0.5)

# Simulate inbred genotypes: 0 or 2, with probability by maf
# P(genotype=2) ~ maf, P(genotype=0) ~ 1-maf; (simple model for inbreds)
G <- matrix(0, nrow = N, ncol = M)
for (j in 1:M) {
  G[, j] <- rbinom(N, size = 1, prob = maf[j]) * 2  # 0 or 2
}
colnames(G) <- paste0("SNP", 1:M)
rownames(G) <- paste0("P", sprintf("%02d", 1:N))

# Quick MAF check (empirical)
emp_maf <- colMeans(G)/2
summary(emp_maf)
```

---

# 2) Simulate true additive marker effects and parent phenotypes

We draw **true additive effects** for markers and compute the **true genetic value** (TGV) for each parent as `G_centered %*% beta` with some noise to create **observed phenotypes**. This becomes the training data for GBLUP.

```r
# Center genotypes to mean 0 for modeling (common for rrBLUP)
G_centered <- scale(G, center = TRUE, scale = FALSE)

# Simulate sparse-ish additive effects
p_causal <- 0.2
beta <- rnorm(M, mean = 0, sd = 0.05) * rbinom(M, 1, p_causal)

# True genetic value (additive)
TGV <- as.vector(G_centered %*% beta)

# Add environmental noise
h2_true <- 0.6
var_g <- var(TGV)
var_e <- var_g * (1 - h2_true)/h2_true
Epsilon <- rnorm(N, mean = 0, sd = sqrt(var_e))

# Observed phenotypes for parents
Y <- TGV + Epsilon

parents <- data.frame(ID = rownames(G), Phenotype = Y, TGV = TGV)
head(parents)
summary(parents$Phenotype)
```

> **Interpretation:** We will try to recover the genetic signal from phenotypes using GBLUP and then predict hybrids from the parental GEBVs.

---

# 3) Train additive GBLUP and compute GEBVs for parents

We fit an **additive GBLUP** using `rrBLUP::mixed.solve` with the marker matrix as a random effect design. It returns **marker effects** `u`, from which we compute **GEBV = G_centered %*% u` for each parent.

```r
# rrBLUP expects Z = marker design (individuals x markers)
fit <- mixed.solve(y = parents$Phenotype, Z = G_centered)

marker_effects <- as.vector(fit$u)            # effects for each SNP
names(marker_effects) <- colnames(G_centered)

GEBV_parents <- as.vector(G_centered %*% marker_effects)
names(GEBV_parents) <- rownames(G_centered)

gebv_df <- data.frame(ID = names(GEBV_parents), GEBV = GEBV_parents) %>%
  left_join(parents, by = "ID")

head(gebv_df)
cor(gebv_df$GEBV, gebv_df$Phenotype)  # sanity check
```

> **Note:** The correlation above gives a rough sense of how well GBLUP recovered the (noisy) phenotypes.

---

# 4) Predict hybrid performance for all possible crosses (additive model)

For inbred parents, a simple **additive-only** hybrid prediction is the **mid-parent value**:
\(
\hat{H}_{ij} = \frac{\text{GEBV}_i + \text{GEBV}_j}{2}
\)

We compute all pairwise crosses \(i<j\).

```r
parents_vec <- gebv_df$ID
K <- length(parents_vec)

# All unique crosses i<j
pairs <- t(combn(parents_vec, 2)) %>% as.data.frame()
colnames(pairs) <- c("Parent1", "Parent2")

# Predicted hybrid value under additive model
id_to_gebv <- setNames(gebv_df$GEBV, gebv_df$ID)
pairs$Hybrid_pred <- (id_to_gebv[pairs$Parent1] + id_to_gebv[pairs$Parent2]) / 2

# Compute parental baselines (using observed phenotypes for heterosis reference)
id_to_pheno <- setNames(gebv_df$Phenotype, gebv_df$ID)
pairs$MP <- (id_to_pheno[pairs$Parent1] + id_to_pheno[pairs$Parent2]) / 2
pairs$BP <- pmax(id_to_pheno[pairs$Parent1], id_to_pheno[pairs$Parent2])

# Heterosis metrics (relative)
pairs$MPH <- (pairs$Hybrid_pred - pairs$MP) / pairs$MP
pairs$BPH <- (pairs$Hybrid_pred - pairs$BP) / pairs$BP

# Also keep absolute versions (optional)
pairs$MPH_abs <- pairs$Hybrid_pred - pairs$MP
pairs$BPH_abs <- pairs$Hybrid_pred - pairs$BP

head(pairs)
summary(pairs$MPH)
```

> **Tip:** If you prefer a fully “predicted” pipeline, you can use **parental GEBVs** instead of observed phenotypes in MP and BP. Here we use observed parental phenotypes as the baseline, which is common in practice.

---

# 5) Visualize the crossing matrix as a heterosis heatmap

We’ll build an \(N \times N\) matrix where cell \((i,j)\) contains the predicted **MPH**. We’ll show only the upper triangle since crosses are symmetric.

```r
# Prepare a square matrix of MPH
mp <- matrix(NA, nrow = N, ncol = N, dimnames = list(parents_vec, parents_vec))
for (k in 1:nrow(pairs)) {
  i <- pairs$Parent1[k]; j <- pairs$Parent2[k]
  mp[i, j] <- pairs$MPH[k]
}

# Melt to long format for ggplot
mp_long <- melt(mp, varnames = c("Parent1", "Parent2"), value.name = "MPH") %>%
  filter(!is.na(MPH))

ggplot(mp_long, aes(x = Parent1, y = Parent2, fill = MPH)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Predicted Mid-Parent Heterosis (MPH) Heatmap",
       x = "Parent 1", y = "Parent 2", fill = "MPH")
```

---

# 6) Rank and display the best crosses

We rank crosses by **MPH** (and also show BPH). You can change `top_k` to any number.

```r
top_k <- 15
best_crosses <- pairs %>%
  arrange(desc(MPH)) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, Parent1, Parent2, Hybrid_pred, MP, BP, MPH, BPH, MPH_abs, BPH_abs) %>%
  head(top_k)

best_crosses
```

---

# 7) (Optional) Add a quick dominance-aware approximation

Dominance effects contribute to heterosis via **heterozygosity** in the hybrid. A simple approximation is:

\(
\hat{H}^{(A+D)}_{ij} \approx \frac{GEBV_i + GEBV_j}{2} + \bar{d} \cdot H_{ij}
\)

where \(H_{ij}\) is the **proportion of loci** where the two parents differ (i.e., expected heterozygote in the F1 if parents are inbred) and \(\bar{d}\) is an **average dominance effect** estimated from data (here we estimate it by regressing `Phenotype` residuals on heterozygosity). This is **very simplified**, but useful for teaching intuition.

```r
# 7.1) Compute pairwise expected heterozygosity between inbred parents
# For inbreds with genotypes 0/2, hybrid is heterozygote when parents differ at a locus.
het_prop <- function(g1, g2) mean(abs(g1 - g2) == 2)

Hprop <- matrix(NA, nrow = N, ncol = N, dimnames = list(parents_vec, parents_vec))
for (i in 1:N) for (j in 1:N) {
  if (i != j) Hprop[i, j] <- het_prop(G[i, ], G[j, ])
}

# 7.2) Estimate an average dominance effect (bar_d) from data (very rough!)
# Regress phenotype residual (after removing additive GEBV) on individual's inbreeding (which is ~0 here)
# Instead, we use cross-level regression: (Hybrid_pred - MP based on GEBV) ~ Hprop
# Build a data frame for all i<j
pairs2 <- t(combn(parents_vec, 2)) %>% as.data.frame()
colnames(pairs2) <- c("Parent1", "Parent2")
pairs2$Hprop <- mapply(function(i,j) Hprop[i, j], pairs2$Parent1, pairs2$Parent2)

# Use additive-only hybrid based on GEBV for both hybrid and MP baselines
pairs2$Hybrid_A <- (id_to_gebv[pairs2$Parent1] + id_to_gebv[pairs2$Parent2]) / 2
pairs2$MP_A     <- (id_to_gebv[pairs2$Parent1] + id_to_gebv[pairs2$Parent2]) / 2

# Empirically regress (Hybrid_A - MP_A) on Hprop: the LHS is 0 by construction,
# so we approximate bar_d simply from the relationship with observed phenotypes as a proxy.
# For a toy illustration, we estimate bar_d by relating (observed mid-parent phenotype) residuals to Hprop.
pairs2$MP_obs <- (id_to_pheno[pairs2$Parent1] + id_to_pheno[pairs2$Parent2]) / 2
lm_d <- lm((pairs2$MP_obs - mean(parents$Phenotype)) ~ pairs2$Hprop)
bar_d <- coef(lm_d)[2]
bar_d
```

With this (very crude) \(\bar{d}\), we can create a **dominance-adjusted** hybrid prediction and recompute heterosis.

```r
pairs2$Hybrid_AD <- pairs2$Hybrid_A + bar_d * pairs2$Hprop

# Recompute heterosis relative to observed parents
pairs2$BP_obs <- pmax(id_to_pheno[pairs2$Parent1], id_to_pheno[pairs2$Parent2])
pairs2$MPH_AD <- (pairs2$Hybrid_AD - pairs2$MP_obs) / pairs2$MP_obs
pairs2$BPH_AD <- (pairs2$Hybrid_AD - pairs2$BP_obs) / pairs2$BP_obs

best_crosses_AD <- pairs2 %>%
  arrange(desc(MPH_AD)) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, Parent1, Parent2, Hybrid_AD, MP_obs, BP_obs, MPH_AD, BPH_AD) %>%
  head(10)

best_crosses_AD
```

> ⚠️ **Caveat:** This dominance section is **for intuition only**. For serious work, use explicit additive + dominance models (e.g., `sommer::mmer` with additive and dominance kernels) or hybrid performance models with dominance terms at the marker level.

---

# 8) What to hand in (for students)

1. Correlation between **GEBV** and **phenotypes** across parents.  
2. **Heatmap** of predicted **MPH** for all crosses.  
3. A table of the **Top 10 crosses** by **MPH** (and BPH if desired).  
4. (Optional) Repeat using the **dominance-aware** approximation and compare rankings.

---

# 9) Extensions & discussion prompts

- Cross-validation: split parents into train/test; evaluate prediction of test phenotypes by GBLUP.  
- Try different `M` and `N` to see how marker density and sample size affect accuracy.  
- Replace observed parental phenotypes in MP/BP with **parental GEBVs** to see differences.  
- Discuss when **BPH** is a more relevant metric than **MPH**.

---

# Session info

```r
sessionInfo()
```
