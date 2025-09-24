# Extra Exercise: Genomic Selection with Random Forest

In this exercise, we will use **genomic selection (GS)** to predict sucrose content from genome-wide SNP markers.  
Unlike the PCA classification exercise (High vs Low sucrose), this is a **regression task**: we predict the actual sucrose values.

---
## Step 0 - Installing the `randomForest` package in RStudio**  

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

---

## **Step 1 — Input data**

We start from:
- **Genotype data**: PLINK `.raw` file (`geno_matrix.raw`), generated from QC’d SNPs.  
- **Phenotype data**: CSV file with sucrose values (`gwas_phen_table_120.csv`).  

```r
library(randomForest)

# Load genotype matrix
geno <- read.table("/lisc/scratch/course/pgbiow/06_diversity_structure/plink/geno_matrix.raw",
                   header=TRUE, check.names=FALSE)

# Load phenotypes (IID + SUC)
pheno <- read.csv("/lisc/scratch/course/pgbiow/data/metadata/gwas_phen_table_120.csv")

# Merge by IID
d <- merge(geno, pheno, by="IID")
d <- d[!is.na(d$SUC), ]

cat("Samples in genotype: ", nrow(geno),  "\n")
cat("Samples in phenotype:", nrow(pheno), "\n")
cat("Samples after merge:", nrow(d), "\n")
```

---

## **Step 2 — SNP filtering and imputation**

PLINK `.raw` files may contain missing values. Random Forest cannot handle `NA`s, so we will filter and impute.  

```r
# Identify SNP columns (skip first 6: FID, IID, PAT, MAT, SEX, PHENOTYPE)
snp_names <- colnames(geno)[7:ncol(geno)]
snp_names <- intersect(snp_names, colnames(d))  # safety check

X <- d[, snp_names, drop=FALSE]

# Filter SNPs with <= 10% missing data
miss_rate <- colMeans(is.na(X))
keep <- miss_rate <= 0.10
X <- X[, keep, drop=FALSE]
cat("SNPs after missingness filter:", ncol(X), "\n")

# Impute remaining missing values with SNP mean
for (j in seq_len(ncol(X))) {
  if (anyNA(X[[j]])) {
    m <- mean(X[[j]], na.rm=TRUE)
    X[[j]][is.na(X[[j]])] <- m
  }
}

# Remove near-constant SNPs (zero variance)
vars <- sapply(X, var)
X <- X[, vars > 0, drop=FALSE]
cat("SNPs after variance filter:", ncol(X), "\n")
```

---

## **Step 3 — Train/test split**

We split the dataset into **70% training** and **30% testing**.  

```r
Y <- d$SUC
IDs <- d$IID

set.seed(42)
train_idx <- sample(seq_len(nrow(X)), size=floor(0.7*nrow(X)))
X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
y_train <- Y[train_idx]
y_test  <- Y[-train_idx]
id_test <- IDs[-train_idx]

cat("Training samples:", nrow(X_train), "\n")
cat("Testing samples:", nrow(X_test), "\n")
```

---

## **Step 4 — Train Random Forest model**

```r
rf <- randomForest(x=X_train, y=y_train,
                   ntree=500, importance=TRUE)
```

---

## **Step 5 — Predictions and accuracy**

```r
pred <- predict(rf, X_test)

# Evaluate accuracy
corr <- cor(pred, y_test, use="complete.obs")
r2   <- corr^2
mae  <- mean(abs(pred - y_test))
rmse <- sqrt(mean((pred - y_test)^2))

cat("Prediction correlation (r):", round(corr, 3), "\n")
cat("Prediction R^2:           ", round(r2,   3), "\n")
cat("MAE:                       ", round(mae,  3), "\n")
cat("RMSE:                      ", round(rmse, 3), "\n")
```

---

## **Step 6 — Save outputs**

```r
# Save predictions
out <- data.frame(IID=id_test, Observed=y_test, Predicted=pred)
write.csv(out, "/lisc/scratch/course/pgbiow/06_diversity_structure/plink/rf_predictions.csv",
          row.names=FALSE)

# Variable importance
vip <- importance(rf)
write.csv(vip, "/lisc/scratch/course/pgbiow/06_diversity_structure/plink/rf_variable_importance.csv")

# Importance plot
png("/lisc/scratch/course/pgbiow/06_diversity_structure/plink/rf_varImp.png",
    width=1000, height=800, res=150)
varImpPlot(rf, main="Random Forest Variable Importance (SNPs)")
dev.off()
```

---

## **Step 7 — Interpretation**

- **Correlation (r)**: how well predictions track observed values.  
- **R²**: proportion of variance in sucrose explained by SNPs.  
- **MAE / RMSE**: prediction error (average absolute vs squared error).  

For selection: the **ranking of predicted values** is what matters.  
Individuals with the **highest predicted sucrose** are the best candidates.  

Example:
```r
top10 <- out[order(-out$Predicted), ][1:10, ]
print(top10)
```

This shows the **top 10 candidates** for sucrose selection based on genomic prediction.
