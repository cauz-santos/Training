## Week 2: Applied Bioinformatics for Genomics and Breeding

## Day 7: SNP Applications in Breeding - GWAS

Day 7 focuses on the practical applications of SNPs in agricultural and animal breeding. We will explore how genetic markers, particularly SNPs, are used to understand population structure, identify desirable traits, and accelerate breeding programs.

**Morning Session 1: Introduction to Breeding Needs**

This session will provide an overview of the challenges and objectives in modern breeding programs and how genomic tools, particularly SNP markers, have revolutionized the field. Traditional breeding relies on phenotypic selection, which can be slow, labor-intensive, and influenced by environmental factors. Genomic breeding aims to accelerate this process by using genetic information directly.

**Key Breeding Objectives:**

*   **Increased Yield:** Developing varieties or breeds that produce more food, fiber, or other products.
*   **Improved Quality:** Enhancing nutritional content, taste, shelf-life, or other desirable product characteristics.
*   **Disease and Pest Resistance:** Breeding for resilience against common pathogens and pests to reduce crop losses and reliance on chemical treatments.
*   **Stress Tolerance:** Developing varieties that can withstand adverse environmental conditions like drought, salinity, or extreme temperatures.
*   **Adaptation to Climate Change:** Creating resilient crops and livestock that can thrive in changing climatic conditions.
*   **Reduced Environmental Impact:** Breeding for traits that require fewer inputs (e.g., water, fertilizer) or produce less waste.

**How Genomics Addresses Breeding Needs:**

Genomic technologies, especially high-throughput SNP genotyping, provide a powerful means to:

*   **Accelerate Selection:** Identify individuals with desirable genes at an early stage (e.g., seedling stage in plants), significantly shortening breeding cycles.
*   **Increase Precision:** Select for complex traits influenced by many genes (quantitative traits) with greater accuracy.
*   **Broaden Genetic Diversity:** Efficiently manage and utilize genetic diversity from germplasm collections.
*   **Understand Genetic Architecture:** Map genes responsible for specific traits, leading to a deeper understanding of the underlying biology.
*   **Overcome Environmental Variation:** Genomic selection is less affected by environmental fluctuations compared to phenotypic selection.

**Morning Session 2: Population Structure and PCA (e.g., PLINK, SNPRelate)**

Understanding the genetic structure of a population is crucial in breeding and genetic studies. Population structure refers to the presence of subgroups within a larger population that have different allele frequencies. Ignoring population structure can lead to spurious associations in genetic studies (false positives) and inefficient breeding strategies.

**What is Population Structure?**

Population structure arises from factors like geographical isolation, limited gene flow, historical migrations, and selection pressures. It can manifest as distinct genetic clusters within a population.

**Principal Component Analysis (PCA) for Population Structure:**

PCA is a widely used statistical method to reduce the dimensionality of complex datasets while retaining most of the variation. In genomics, PCA is applied to SNP data to visualize genetic relationships among individuals and identify population structure. Individuals that are genetically similar will cluster together in a PCA plot.

*   **How it works:** PCA transforms the original SNP data into a new set of uncorrelated variables called principal components (PCs). The first few PCs capture the largest proportion of genetic variation in the dataset.
*   **Interpretation:** By plotting individuals based on their scores on the first two or three PCs, we can often observe distinct clusters corresponding to different populations or ancestral groups.

**Tools for Population Structure Analysis:**

    *   **PLINK:** A comprehensive open-source whole-genome association analysis toolset. It can perform various genetic data manipulations, including calculating principal components for population structure analysis.
        *   **Input:** Typically requires genotype data in PLINK binary format (BED, BIM, FAM).
        *   **Command for PCA:**
            ```bash
            # Assuming you have a PLINK binary file set (my_data.bed, my_data.bim, my_data.fam)
            plink --bfile my_data --pca --out my_pca_results
            # Output: my_pca_results.eigenvec (eigenvectors/PC scores) and my_pca_results.eigenval (eigenvalues)
            ```
        *   **Output:** Files containing eigenvector (PC scores for each individual) and eigenvalue (proportion of variance explained by each PC) information.
    *   **SNPRelate (R package):** An R package for genome-wide association studies (GWAS) and population genetics. It provides functions for calculating genetic relatedness, principal component analysis, and other population genetic analyses.
        *   **Input:** GDS (Genomic Data Structure) format, which can be created from VCF files.
        *   **Example R code:**
            ```R
            library(SNPRelate)

            # Convert VCF to GDS format (if not already in GDS)
            # snpgdsVCF2GDS("input.vcf", "output.gds", method="biallelic.only")

            # Open the GDS file
            genofile <- snpgdsOpen("output.gds")

            # Perform PCA
            pca <- snpgdsPCA(genofile, num.thread=2)

            # Get sample IDs and population information (if available)
            sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
            # pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop_code")) # If population info is in GDS

            # Plotting the PCA results
            plot(pca$eigenvect[,1], pca$eigenvect[,2], 
                 xlab=paste0("PC1 (", round(pca$varprop[1]*100, 2), "%)"), 
                 ylab=paste0("PC2 (", round(pca$varprop[2]*100, 2), "%)"),
                 main="PCA of SNP Data")
            # Add text labels or color by population if pop_code is available
            # text(pca$eigenvect[,1], pca$eigenvect[,2], labels=sample.id, cex=0.7, pos=3)

            # Close the GDS file
            snpgdsClose(genofile)
            ```

**Morning Session 3: Introgression and Admixture Analysis (e.g., Treemix, Dsuite)**

Introgression and admixture are evolutionary processes that involve the mixing of genetic material between previously divergent populations or species. Understanding these processes is vital for studying evolutionary history, identifying adaptive traits, and managing genetic resources in breeding.

**Introgression:** The movement of genes from one species or population into another through hybridization and backcrossing.

**Admixture:** The result of interbreeding between two or more previously isolated populations, leading to individuals with mixed ancestry.

**Tools for Introgression and Admixture Analysis:**

    *   **TreeMix:** A powerful program for inferring population trees that include migration events (gene flow/admixture). It uses allele frequency data to build a maximum likelihood tree and then adds migration edges to improve the fit of the tree to the data.
        *   **Input:** Allele frequency data (e.g., from VCF files, converted to TreeMix format using `vcf2treemix.py` or similar).
        *   **Output:** Newick tree file with migration edges, and residual plots.
        *   **Example Workflow (Conceptual):**
            ```bash
            # 1. Convert VCF to TreeMix input format (requires a script like vcf2treemix.py)
            # python vcf2treemix.py -i your_variants.vcf.gz -o treemix_input.gz

            # 2. Run TreeMix (allowing for 3 migration events)
            treemix -i treemix_input.gz -m 3 -o output_treemix

            # 3. Plot the tree (requires R and the TreeMix plotting script)
            # Rscript plot_tree.r output_treemix
            ```
    *   **Dsuite:** A suite of tools for calculating D-statistics (also known as ABBA-BABA test) and related statistics to detect introgression. D-statistics measure the excess of shared derived alleles between non-sister lineages, indicating gene flow.
        *   **Input:** VCF file.
        *   **Output:** Tables with D-statistics and Z-scores.
        *   **Example Workflow (Conceptual):**
            ```bash
            # 1. Create a population file (e.g., populations.txt) with SampleID\tPopulationName
            # 2. Run Dsuite Dtrios to calculate D-statistics for all trios of populations
            Dsuite Dtrios my_data.vcf.gz populations.txt

            # 3. Run Dsuite Dinvestigate to test for introgression from a specific outgroup
            # Dsuite Dinvestigate my_data.vcf.gz populations.txt outgroup_population
            ```

**Afternoon Session 4: GWAS and Genotype-Phenotype Association Studies (e.g., PLINK)**

Genome-Wide Association Studies (GWAS) are a powerful approach to identify genetic variants (typically SNPs) that are statistically associated with a particular trait or disease. GWAS has been widely applied in human genetics, agriculture, and animal breeding to discover genetic markers linked to complex traits.

**The Principle of GWAS:**

GWAS involves scanning the entire genome for common genetic variants that occur more frequently in individuals with a particular trait (e.g., high yield, disease resistance) compared to individuals without that trait. It is based on the concept of linkage disequilibrium (LD), where alleles at different loci are inherited together more often than expected by chance.

**Steps in a GWAS:**

1.  **Phenotype Collection:** Accurate and precise measurement of the trait of interest across a large number of individuals.
2.  **Genotype Data:** High-throughput genotyping of SNPs across the genome for all individuals.
3.  **Quality Control:** Rigorous QC of both genotype and phenotype data to remove errors and biases.
4.  **Statistical Association:** Performing statistical tests (e.g., logistic regression for binary traits, linear regression for quantitative traits) to assess the association between each SNP and the trait.
5.  **Correction for Multiple Testing:** Adjusting p-values for the large number of statistical tests performed across the genome (e.g., Bonferroni correction, False Discovery Rate).
6.  **Interpretation and Follow-up:** Identifying significant SNPs and exploring their biological relevance.

**Tools for GWAS:**

    *   **PLINK:** As mentioned before, PLINK is a versatile tool for GWAS. It can perform various association tests, including basic case-control association, quantitative trait association, and more complex models incorporating covariates.
        *   **Input:** Genotype data (BED, BIM, FAM) and phenotype data (a simple text file with individual IDs and trait values).
        *   **Example Command for Quantitative Trait Association (Linear Regression):**
            ```bash
            # Assuming you have a PLINK binary file set (my_genotypes.bed, .bim, .fam)
            # And a phenotype file (my_phenotypes.txt) with FID, IID, and phenotype column (e.g., TRAIT1)
            # Example my_phenotypes.txt:
            # FID IID TRAIT1
            # Sample1 Sample1 10.5
            # Sample2 Sample2 12.3

            plink --bfile my_genotypes --pheno my_phenotypes.txt --pheno-name TRAIT1 --linear --out my_gwas_results_linear
            # --linear: Performs a linear regression association test.
            # --pheno-name: Specifies the phenotype column name if multiple are present.
            # Output: my_gwas_results_linear.assoc.linear (contains association results)
            ```
        *   **Example Command for Case-Control Association (Logistic Regression):**
            ```bash
            # Assuming phenotype file has 1 for controls, 2 for cases (or -9 for missing)
            plink --bfile my_genotypes --pheno my_phenotypes.txt --pheno-name DISEASE --logistic --out my_gwas_results_logistic
            # --logistic: Performs a logistic regression association test.
            # Output: my_gwas_results_logistic.assoc.logistic (contains association results)
            ```
        *   **Output:** Results files containing SNP information, allele frequencies, p-values, and effect sizes.
    *   **GEMMA, GCTA, TASSEL:** Other popular tools for GWAS, often offering more advanced statistical models or specific functionalities.
        *   **GEMMA (Genome-wide Efficient Mixed Model Association):** Used for mixed model association mapping, which accounts for population structure and relatedness.
            ```bash
            # Example GEMMA workflow (conceptual)
            # 1. Prepare input files (binary PLINK format and phenotype file)
            # 2. Calculate kinship matrix
            # gemma -bfile my_genotypes -gk 1 -o my_kinship

            # 3. Run association analysis with mixed model
            # gemma -bfile my_genotypes -k output/my_kinship.cXX.txt -lmm 1 -o my_gemma_results
            ```

**Afternoon Session 5: Marker-Assisted Selection, Genomic Selection, and Prediction**

Building upon GWAS, this session will introduce advanced breeding strategies that leverage genomic information to make more efficient and accurate selection decisions.

**Marker-Assisted Selection (MAS):**

MAS involves using DNA markers (like SNPs) that are tightly linked to genes controlling desirable traits to select individuals in a breeding program. Instead of waiting for the phenotype to express, breeders can select based on the presence of the marker allele.

*   **Advantages:** Faster, more efficient, less influenced by environment, can select for traits difficult to phenotype.
*   **Limitations:** Requires markers to be in strong linkage disequilibrium with the causal gene, primarily effective for traits controlled by a few major genes.

**Genomic Selection (GS) and Prediction:**

Genomic selection is a more advanced form of MAS, particularly suited for complex traits controlled by many genes with small effects. Instead of using a few major markers, GS uses all available markers across the genome to predict an individual\`s breeding value (Genomic Estimated Breeding Value - GEBV).

*   **Training Population:** A population of individuals that have both genotype and phenotype data. This population is used to estimate the effects of all markers across the genome.
*   **Prediction Equation:** A statistical model is built using the training population to predict GEBVs for new individuals.
*   **Selection Candidate Population:** Individuals that are genotyped but not yet phenotyped (or are too young to be phenotyped). Their GEBVs are predicted using the equation from the training population.
*   **Advantages:** More accurate for complex traits, can accelerate breeding cycles significantly, no need to identify individual causal genes.
*   **Tools:** Software like `GBLUP`, `rrBLUP` (R packages), and `AlphaSimR` (for simulation) are used for genomic prediction.

**Conceptual Workflow of Genomic Selection:**

1.  **Genotype and Phenotype Training Population:** Collect SNP data and trait data for a diverse set of individuals.
2.  **Build Prediction Model:** Use statistical methods to estimate the effect of each SNP on the trait.
3.  **Genotype Selection Candidates:** Obtain SNP data for individuals in the breeding population that need to be selected.
4.  **Predict GEBVs:** Use the established model to predict the breeding value for each selection candidate.
5.  **Select Best Individuals:** Choose individuals with the highest GEBVs for the next breeding cycle.

By the end of Day 7, you will have a comprehensive understanding of how genomic data, particularly SNPs, are applied in modern breeding programs to enhance efficiency, precision, and genetic gain. You will be familiar with key concepts like population structure, introgression, GWAS, MAS, and genomic selection.


