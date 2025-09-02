## Week 2: Applied Bioinformatics for Genomics and Breeding

### Day 9: RNA-seq Data Analysis

Day 9 is dedicated to RNA sequencing (RNA-seq) data analysis, a powerful high-throughput technology used to measure gene expression levels, discover novel transcripts, and identify gene fusions. We will cover the entire workflow, from experimental design and pre-processing to mapping, quantification, and differential expression analysis.

**Morning Session 1: RNA-seq Experimental Design and Pre-processing**

As discussed on Day 5, robust experimental design is critical for RNA-seq studies. This session will reiterate key design principles specific to RNA-seq and then dive into the initial steps of data pre-processing.

**Key Considerations for RNA-seq Experimental Design:**

*   **Biological Replicates:** Absolutely essential for statistical power and distinguishing true biological variation from technical noise. Aim for at least 3-5 biological replicates per condition.
*   **Sequencing Depth:** The number of reads per sample. Deeper sequencing allows for better detection of lowly expressed genes and more accurate quantification. Typical depths range from 10-50 million reads per sample, depending on the complexity of the transcriptome and research question.
*   **Read Length:** Longer reads can improve alignment accuracy, especially in regions with repetitive sequences or alternative splicing.
*   **Paired-end vs. Single-end:** Paired-end reads provide more information for alignment and transcript assembly, especially for identifying splice junctions.
*   **RNA Quality and Integrity:** High-quality RNA is crucial. RNA Integrity Number (RIN) scores are commonly used to assess RNA degradation. Degraded RNA can lead to biased results.
*   **Batch Effects:** As mentioned, samples processed at different times or by different personnel can introduce systematic biases. Randomization of samples across batches is important.

**RNA-seq Pre-processing Steps:**

After sequencing, raw RNA-seq reads (FASTQ files) need to be pre-processed to remove low-quality data and artifacts before alignment and quantification.

1.  **Quality Control (QC):** Use tools like FastQC (covered on Day 5) to assess the quality of raw reads. Look for issues like low quality scores, adapter contamination, and GC content bias.
2.  **Trimming and Filtering:** Use tools like Trimmomatic (covered on Day 5) or Cutadapt to remove adapter sequences, low-quality bases, and short reads. This step is crucial for improving alignment accuracy and reducing computational burden.

    ```bash
    # Example using Cutadapt (alternative to Trimmomatic)
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
             -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
             raw_reads_R1.fastq.gz raw_reads_R2.fastq.gz
    ```

**Morning Session 2: Mapping and Quantification (e.g., STAR)**

Once reads are pre-processed, the next step is to align them to a reference genome and then quantify gene or transcript expression levels.

**Read Mapping (Alignment):**

RNA-seq reads are aligned to a reference genome to determine their genomic origin. Unlike DNA sequencing, RNA-seq reads originate from transcribed regions (exons) and often span exon-exon junctions. Therefore, specialized splice-aware aligners are required.

*   **STAR (Spliced Transcripts Alignment to a Reference):** One of the fastest and most accurate splice-aware aligners for RNA-seq reads. It uses an uncompressed suffix array index to achieve high mapping speeds.

    **STAR Workflow:**

    1.  **Genome Indexing:** Build a STAR genome index from the reference genome FASTA file and gene annotation (GTF/GFF3).
        ```bash
        STAR --runMode genomeGenerate \
             --genomeDir /path/to/STAR_index \
             --genomeFastaFiles /path/to/reference.fasta \
             --sjdbGTFfile /path/to/annotation.gtf \
             --sjdbOverhang 100 # (ReadLength - 1)
        ```
    2.  **Read Alignment:** Align trimmed FASTQ reads to the generated index.
        ```bash
        STAR --runMode alignReads \
             --genomeDir /path/to/STAR_index \
             --readFilesIn trimmed_R1.fastq.gz trimmed_R2.fastq.gz \
             --readFilesCommand zcat \
             --outFileNamePrefix sample_ \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard
        ```

**Quantification:**

After alignment, the next step is to quantify how many reads map to each gene or transcript. This provides a measure of gene expression.

*   **FeatureCounts:** A highly efficient and flexible program for counting reads to genomic features (genes, exons, etc.) from BAM files. It is often used after STAR alignment.

    ```bash
    featureCounts -p -t exon -g gene_id \
                  -a /path/to/annotation.gtf \
                  -o gene_counts.txt \
                  sample_Aligned.sortedByCoord.out.bam
    ```
    *   `-p`: Paired-end reads.
    *   `-t exon`: Count reads mapping to exons.
    *   `-g gene_id`: Group counts by gene ID.
    *   `-a`: Path to the GTF annotation file.

**Morning Session 3: Differential Expression Analysis in R (DESeq2, edgeR)**

The ultimate goal of many RNA-seq experiments is to identify genes that are differentially expressed (DEGs) between different conditions (e.g., treated vs. control, disease vs. healthy). This involves statistical modeling to account for variability and normalization.

**Normalization:**

Raw read counts cannot be directly compared between samples due to differences in sequencing depth, gene length, and RNA composition. Normalization methods adjust raw counts to make them comparable. Common methods include TPM (Transcripts Per Million), FPKM (Fragments Per Kilobase of transcript per Million mapped reads), and methods used by DESeq2/edgeR.

**Differential Expression Analysis (DEA):**

DEA uses statistical models (often based on negative binomial distribution) to identify genes whose expression levels are significantly different between experimental groups.

**R Packages for DEA:**

    *   **DESeq2:** A popular Bioconductor package for differential gene expression analysis based on the negative binomial distribution. It is robust to outliers and provides good control over false discovery rates.

    **DESeq2 Workflow (Conceptual):**

    1.  **Load Data:** Read gene count matrix (e.g., from FeatureCounts) and sample metadata into R.
    2.  **Create DESeqDataSet Object:** Combine counts and metadata.
    3.  **Run DESeq:** Perform normalization and differential expression testing.
    4.  **Extract Results:** Get a table of differentially expressed genes with p-values, adjusted p-values, and fold changes.
    5.  **Visualization:** Create MA plots, volcano plots, and heatmaps to visualize results.

    ```R
    # Example R code for DESeq2
    library(DESeq2)

    # Assuming \'counts_matrix\' is your gene counts and \'colData\' is your sample metadata
    # colData should be a data.frame with row names matching column names of counts_matrix
    # and a column indicating experimental condition (e.g., \'condition\')

    # Create a dummy counts matrix and colData for demonstration
    # In a real scenario, counts_matrix would come from featureCounts output
    # and colData from your experimental design file.
    set.seed(123)
    counts_matrix <- matrix(rnbinom(n=100*6, size=10, mu=100), ncol=6)
    rownames(counts_matrix) <- paste0("gene", 1:100)
    colnames(counts_matrix) <- paste0("sample", 1:6)

    colData <- data.frame(condition = factor(c("control", "control", "control", "treated", "treated", "treated")))
    rownames(colData) <- colnames(counts_matrix)

    dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                  colData = colData,
                                  design = ~ condition)

    # Pre-filtering (optional but recommended for large datasets)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]

    dds <- DESeq(dds)
    res <- results(dds)

    # Order by adjusted p-value
    res_ordered <- res[order(res$padj),]

    # Summary of results
    summary(res)

    # Plotting (e.g., MA plot)
    plotMA(res, main="MA Plot")

    # Export results
    # write.csv(as.data.frame(res_ordered), file="deseq2_results.csv")
    ```

*   **edgeR:** Another widely used Bioconductor package for differential expression analysis of RNA-seq data, also based on the negative binomial distribution. It offers similar functionalities to DESeq2.

    **edgeR Workflow (Conceptual):**

    1.  **Load Data:** Read gene count matrix and sample metadata into R.
    2.  **Create DGEList Object:** Create an edgeR DGEList object from counts and group information.
    3.  **Normalization:** Perform TMM (Trimmed Mean of M-values) normalization.
    4.  **Estimate Dispersion:** Estimate common, trended, and tagwise dispersions.
    5.  **Fit Model and Test:** Fit a negative binomial generalized linear model and perform statistical tests for differential expression.
    6.  **Extract Results:** Get a table of differentially expressed genes.

    ```R
    # Example R code for edgeR
    library(edgeR)

    # Assuming \'counts_matrix\' and \'colData\' from DESeq2 example
    group <- colData$condition

    y <- DGEList(counts=counts_matrix, group=group)

    # Pre-filtering
    keep <- filterByExpr(y)
    y <- y[keep,,keep.lib.sizes=FALSE]

    # Normalization
    y <- calcNormFactors(y)

    # Estimate dispersion
    y <- estimateDisp(y)

    # Fit GLM and test for differential expression
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef=2) # Assuming \'conditiontreated\' is the second coefficient

    # Extract results
    results_edgeR <- topTags(qlf, n=nrow(y))
    print(results_edgeR)

    # Plotting (e.g., MDS plot, MA plot)
    # plotMDS(y)
    # plotSmear(qlf, de.tags=rownames(results_edgeR$table))
    ```

**Afternoon Session 4: Discussion: Issues Arising from the Course**

This final session of the course is dedicated to an open discussion, allowing participants to reflect on the topics covered, ask questions, and discuss challenges or issues they encountered during the two weeks. It\'s an opportunity to consolidate learning and address any remaining uncertainties.

**Potential Discussion Topics:**

*   **Troubleshooting:** Common errors encountered during command-line operations, scripting, or tool execution.
*   **Data Interpretation:** Challenges in interpreting QC reports, alignment statistics, variant calls, or differential expression results.
*   **Choosing the Right Tool:** When to use which tool for a specific task (e.g., GATK vs. bcftools for SNP calling, DESeq2 vs. edgeR for DEA).
*   **Computational Resources:** Strategies for managing large datasets and utilizing HPC resources effectively.
*   **Experimental Design Revisited:** How to design future experiments to minimize common pitfalls.
*   **Next Steps:** Resources for continued learning, advanced topics not covered in the course, and how to apply these skills to individual research projects.
*   **Feedback:** Constructive feedback on the course content, pace, and delivery.

This session aims to ensure that all participants leave the course with a clear understanding of the foundational concepts and practical skills, and feel confident in applying them to their own bioinformatics challenges. It also provides an opportunity for networking and building a community of practice among the participants.

**Conclusion:**

This two-week intensive bioinformatics training course has provided a comprehensive overview of essential concepts and practical skills in handling and analyzing genomic data. From mastering the Linux command line to performing advanced RNA-seq analysis, you have gained valuable experience that will serve as a strong foundation for your future endeavors in bioinformatics. I encourage you to continue exploring, practicing, and applying these skills to unlock new biological insights. The field of bioinformatics is constantly evolving, and continuous learning is key to staying at the forefront of discovery. I wish you all the best in your bioinformatics journey!

# Congratulations!!! You have completed the training



## References

[1] The Linux Documentation Project. Available at: https://www.tldp.org/
[2] Python Documentation. Available at: https://docs.python.org/3/
[3] Biopython Project. Available at: https://biopython.org/
[4] R Project for Statistical Computing. Available at: https://www.r-project.org/
[5] NCBI (National Center for Biotechnology Information). Available at: https://www.ncbi.nlm.nih.gov/
[6] European Bioinformatics Institute (EMBL-EBI). Available at: https://www.ebi.ac.uk/
[7] FastQC. Available at: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[8] MultiQC. Available at: https://multiqc.info/
[9] Trimmomatic. Available at: http://www.usadellab.org/cms/?page=trimmomatic
[10] GATK (Genome Analysis Toolkit). Available at: https://gatk.broadinstitute.org/hc/en-us
[11] Stacks. Available at: http://catchenlab.com/stacks/
[12] bcftools. Available at: http://samtools.github.io/bcftools/
[13] VCFtools. Available at: https://vcftools.github.io/index.html
[14] SnpEff. Available at: http://pcingola.github.io/SnpEff/
[15] PLINK. Available at: https://www.cog-genomics.org/plink/
[16] SNPRelate. Available at: https://bioconductor.org/packages/release/bioc/html/SNPRelate.html
[17] TreeMix. Available at: https://bitbucket.org/nygcresearch/treemix/wiki/Home
[18] Dsuite. Available at: https://github.com/millanek/Dsuite
[19] STAR (Spliced Transcripts Alignment to a Reference). Available at: https://github.com/alexdobin/STAR
[20] FeatureCounts. Available at: http://subread.sourceforge.net/
[21] DESeq2. Available at: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
[22] edgeR. Available at: https://bioconductor.org/packages/release/bioc/html/edgeR.html
[23] scikit-learn. Available at: https://scikit-learn.org/stable/
[24] pandas. Available at: https://pandas.pydata.org/
[25] numpy. Available at: https://numpy.org/
[26] matplotlib. Available at: https://matplotlib.org/
[27] seaborn. Available at: https://seaborn.pydata.org/


