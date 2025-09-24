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


### Step 2 - Input Data

We will use:

- **Parent genotypes** from: `training_sets/parents/parents.raw`
- **Offspring genotypes** from: `training_sets/offsprings/offsprings.raw`
- **Offspring phenotypes** from: `testing_sets/list_input_parent_pairs_offsprings_*.csv`
- **Parent phenotypes** from: `training_sets/ind_names_phenotypes_heterosis_Sf_30.csv`

```{r input-files}
offspring_raw <- file.path(base, "training_sets/offsprings/offsprings.raw")
parents_raw   <- file.path(base, "training_sets/parents/parents.raw")
pairs_csv     <- list.files(file.path(base, "testing_sets"),
                            pattern="^list_input_parent_pairs_offsprings_.*\\.csv$", full.names=TRUE)[1]
parents_pheno_csv <- file.path(base, "training_sets", "ind_names_phenotypes_heterosis_Sf_30.csv")

TRAIT <- "MBW-MA"   # <-- change this to target other traits if needed
```


### Step 3 - Load and Harmonize Data

```{r load-data}
read_raw <- function(path) {
  x <- read_table(path, show_col_types = FALSE, na = c("NA"))
  x %>% mutate(ID = as.character(IID)) %>% select(-FID, -IID, -PAT, -MAT, -SEX, -PHENOTYPE)
}

offs_raw <- read_raw(offspring_raw)
pars_raw <- read_raw(parents_raw)

# Intersect SNPs
snps <- intersect(setdiff(names(offs_raw), "ID"), setdiff(names(pars_raw), "ID"))
offs <- offs_raw %>% select(ID, all_of(snps))
pars <- pars_raw %>% select(ID, all_of(snps))

# Offspring phenotypes
pairs_raw <- read_csv(pairs_csv, show_col_types = FALSE)
offs_pheno <- pairs_raw %>%
  transmute(ID = as.character(id),
            Trait = as.numeric(.data[[TRAIT]]),
            Female = as.character(`Female.Parents`),
            Male   = as.character(`Male.Parents`)) %>%
  filter(!is.na(ID), !is.na(Trait))

offs <- offs %>% semi_join(offs_pheno, by = "ID")
offs_pheno <- offs_pheno %>% semi_join(offs, by = "ID")
```


### Step 4 - Build Additive and Dominance Codes

```{r add-dom-codes}
G_offs <- as.matrix(offs %>% select(-ID))
col_means <- colMeans(G_offs, na.rm = TRUE)
p <- col_means / 2
for (j in seq_len(ncol(G_offs))) {
  G_offs[is.na(G_offs[,j]), j] <- col_means[j]
}
Z_A <- sweep(G_offs, 2, 2*p, "-")
Z_D <- (G_offs == 1) * 1.0
EH  <- 2 * p * (1 - p)
Z_D <- sweep(Z_D, 2, EH, "-")
rownames(Z_A) <- offs$ID; rownames(Z_D) <- offs$ID
```


### Step 5 - Train rrBLUP (Additive + Dominance)

```{r rrblup-train}
y <- offs_pheno$Trait
idx <- match(offs_pheno$ID, rownames(Z_A))
ZA <- Z_A[idx, , drop=FALSE]; ZD <- Z_D[idx, , drop=FALSE]
Z  <- cbind(ZA, ZD)

fit <- mixed.solve(y = y, Z = Z)
beta <- as.vector(fit$u)
beta_A <- beta[seq_len(ncol(ZA))]
beta_D <- beta[seq_len(ncol(ZD)) + ncol(ZA)]
names(beta_A) <- colnames(ZA); names(beta_D) <- colnames(ZD)
```


### Step 6 - Predict Parent × Parent Crosses

```{r predict-crosses}
# Helper function for dominance probabilities
het_prob <- function(g1, g2) {
  if (is.na(g1) || is.na(g2)) return(NA_real_)
  if (g1==g2 && (g1==0 || g1==2)) return(0)
  if ((g1==0 && g2==2) || (g1==2 && g2==0)) return(1)
  if ((g1==1 && g2 %in% c(0,2)) || (g2==1 && g1 %in% c(0,2)) || (g1==1 && g2==1)) return(0.5)
  return(NA_real_)
}
het_prob_vec <- function(v1, v2) vapply(seq_along(v1), function(i) het_prob(v1[i], v2[i]), numeric(1))

G_par <- as.matrix(pars %>% select(-ID))
for (j in seq_len(ncol(G_par))) {
  G_par[is.na(G_par[,j]), j] <- col_means[j]
}
rownames(G_par) <- pars$ID

parents_vec <- rownames(G_par)
pair_idx <- t(combn(seq_along(parents_vec), 2))
Parent1 <- parents_vec[pair_idx[,1]]
Parent2 <- parents_vec[pair_idx[,2]]

A_ij <- (G_par[Parent1, , drop=FALSE] + G_par[Parent2, , drop=FALSE]) / 2
D_ij <- matrix(NA_real_, nrow = nrow(A_ij), ncol = ncol(A_ij),
               dimnames = list(NULL, colnames(A_ij)))
for (k in seq_len(nrow(A_ij))) {
  g1 <- G_par[Parent1[k], ]
  g2 <- G_par[Parent2[k], ]
  D_ij[k, ] <- het_prob_vec(g1, g2)
}

A_ij_c <- sweep(A_ij, 2, 2*p, "-")
D_ij_c <- sweep(D_ij, 2, EH, "-")

hyb_pred <- as.vector(A_ij_c %*% beta_A + D_ij_c %*% beta_D)

pred_df <- tibble(
  Parent1 = Parent1,
  Parent2 = Parent2,
  Hybrid_pred = hyb_pred
)
```

### Step 7 - Compute Heterosis

```{r heterosis}
parents_ph <- read_csv(parents_pheno_csv, show_col_types = FALSE)
id_col <- intersect(names(parents_ph), c("ID","Id","id","name","Name","ind","individual","genotype"))[1]
trait_par <- setdiff(names(parents_ph), c(id_col, grep("heterosis|MPH|BPH|pop|group|population", names(parents_ph), ignore.case=TRUE, value=TRUE)))[1]

parents_ph <- parents_ph %>%
  transmute(ID = as.character(.data[[id_col]]),
            Trait = as.numeric(.data[[trait_par]])) %>%
  distinct()

id2p <- setNames(parents_ph$Trait, parents_ph$ID)

pred_df <- pred_df %>%
  mutate(P1 = id2p[Parent1],
         P2 = id2p[Parent2]) %>%
  filter(!is.na(P1), !is.na(P2)) %>%
  mutate(MP = (P1 + P2)/2,
         BP = pmax(P1, P2),
         MPH = (Hybrid_pred - MP)/MP,
         BPH = (Hybrid_pred - BP)/BP,
         MPH_abs = Hybrid_pred - MP,
         BPH_abs = Hybrid_pred - BP)

best <- pred_df %>% arrange(desc(MPH)) %>%
  mutate(Rank = row_number())
```


### Step 8 - Save Results

```{r save-results}
write_csv(best, file.path(out_dir, paste0("Verdant_crosses_ranked_", TRAIT, ".csv")))

# Heatmap
parents_all <- sort(unique(c(best$Parent1, best$Parent2)))
mat <- matrix(NA_real_, nrow=length(parents_all), ncol=length(parents_all),
              dimnames=list(parents_all, parents_all))
for (k in seq_len(nrow(best))) {
  i <- best$Parent1[k]; j <- best$Parent2[k]
  mat[i, j] <- best$MPH[k]
}
dfh <- as.data.frame(as.table(mat)) %>%
  rename(Parent1=Var1, Parent2=Var2, MPH=Freq) %>%
  filter(!is.na(MPH))

p <- ggplot(dfh, aes(x=Parent1, y=Parent2, fill=MPH)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_minimal(base_size=10) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(title=paste0("Predicted MPH (", TRAIT, ")"), fill="MPH")

ggsave(file.path(out_dir, paste0("MPH_heatmap_", TRAIT, ".png")), p, width=12, height=10, dpi=220)
```


## Conclusion

You have now:

- Trained a **GS model with additive + dominance effects** on the offspring data  
- Predicted **all possible parent × parent crosses**  
- Computed **heterosis measures (MPH, BPH)**  
- Produced a **ranked CSV** and a **heatmap** for visualization  
