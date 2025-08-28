## Week 1: Foundations of Bioinformatics and Genomic Data Handling

## Day 4: Experimental Design, Sequencing and Data QC 

Day 4 is crucial for understanding the importance of data quality in bioinformatics and how experimental design impacts downstream analysis. We will also have a challenging hands-on session to apply our Linux and HPC knowledge.

### Integrative Discussion on Experimental Design and Its Impact on SNP Analysis and RNA-seq

Before diving into data analysis, it is paramount to understand the principles of good experimental design. A well-designed experiment ensures that the data collected is robust, reliable, and capable of answering the biological questions posed. Poor experimental design can lead to biased results, false conclusions, and wasted resources.

**Key Principles of Experimental Design in Genomics:**

*   **Replication:** Having multiple biological replicates is essential to account for biological variability and to ensure statistical power. Without sufficient replicates, it\\'s impossible to distinguish true biological effects from random noise.
*   **Randomization:** Randomly assigning samples to experimental groups or processing order helps to minimize systematic bias.
*   **Blinding:** If possible, blinding researchers to the experimental conditions can prevent unconscious bias during data collection or analysis.
*   **Controls:** Including appropriate positive and negative controls is vital to validate experimental procedures and interpret results correctly.
*   **Sample Size Calculation:** Determining the appropriate number of samples needed to detect a statistically significant effect, considering factors like effect size, variability, and desired statistical power.

**Impact on SNP Analysis:**

For SNP (Single Nucleotide Polymorphism) analysis, experimental design considerations include:

*   **Population Structure:** Understanding the genetic relationships within your sample population is critical to avoid spurious associations. Techniques like Principal Component Analysis (PCA) can help identify population structure.
*   **Relatedness:** Accounting for relatedness among individuals in a study is important to prevent false positives in association studies.
*   **Coverage:** Sufficient sequencing depth (coverage) is needed to accurately call SNPs and distinguish true variants from sequencing errors.
*   **Allele Frequency:** Rare variants require larger sample sizes to achieve statistical significance.

**Impact on RNA-seq:**

For RNA-seq (RNA sequencing), experimental design considerations include:

*   **Biological vs. Technical Replicates:** Biological replicates (different biological samples under the same condition) are essential for statistical inference. Technical replicates (re-sequencing the same sample) are less critical with modern sequencing technologies but can help assess technical variability.
*   **Sequencing Depth:** Adequate sequencing depth is required to accurately quantify gene expression, especially for lowly expressed genes.
*   **Batch Effects:** Samples processed at different times or by different personnel can introduce systematic biases (batch effects). Randomizing samples across batches or accounting for batch effects during analysis is crucial.
*   **RNA Quality:** The quality of input RNA significantly impacts RNA-seq results. RNA integrity (e.g., RIN score) should be assessed.

**Discussion Points:**

*   Share examples of good and bad experimental designs you have encountered.
*   Discuss the challenges of experimental design in real-world biological research.
*   How can bioinformatics tools help identify issues related to experimental design post-data collection?

---

### Understanding Sequencing Technologies and Read Types

Modern genomics relies on different sequencing platforms, each producing data with unique characteristics.  
Choosing the right sequencing technology and read type is critical because it directly impacts the quality of downstream analyses such as **SNP discovery, genome assembly, or transcriptomics**.  

Some technologies produce **short but highly accurate reads** (e.g., Illumina), while others generate **very long reads** that can span complex regions (e.g., PacBio, Oxford Nanopore).  
Similarly, read types such as **single-end, paired-end, long, and ultra-long** each provide different levels of information about the genome or transcriptome.  

In this section, we summarize the most common **sequencing technologies** and **types of reads**, highlighting their advantages, limitations, and typical use cases in bioinformatics and plant breeding.

**Table. Types of Sequencing and Reads**

| Sequencing Type        | Description                                                                 | Typical Applications                          | Read Types & Lengths                |
|------------------------|-----------------------------------------------------------------------------|-----------------------------------------------|--------------------------------------|
| **Sanger Sequencing**  | First-generation sequencing, low throughput but very accurate.              | Small gene fragments, validation, cloning.    | Single reads, ~500–1000 bp           |
| **Illumina (NGS)**     | Short-read sequencing with high accuracy and throughput.                    | WGS, RNA-seq, SNP genotyping, metagenomics.   | Paired-end (2 × 50–300 bp) or single-end |
| **Ion Torrent**        | Short-read sequencing, detects pH changes during DNA synthesis.             | Small genomes, targeted sequencing.           | Single or paired-end, up to ~400 bp  |
| **PacBio (HiFi/CLR)**  | Long-read sequencing; HiFi reads combine length with high accuracy.         | Genome assembly, isoform sequencing, SVs.     | Continuous Long Reads (up to 50 kb) or HiFi (~10–25 kb, >99% accurate) |
| **Oxford Nanopore (ONT)** | Portable, real-time long-read sequencing.                               | Structural variants, genome assembly, direct RNA sequencing. | Reads from kb to Mb (very long), lower raw accuracy but improving with new chemistries |
| **10x Genomics (Linked-Reads)** | Uses barcoding to connect short Illumina reads to long molecules. | Phasing, haplotyping, structural variation.   | Short reads (~150 bp) but linked across long DNA molecules (~50–100 kb) |


**Table. Types of Reads:**

| Read Type        | Description                                                                 | Example Use Case                             |
|------------------|-----------------------------------------------------------------------------|---------------------------------------------|
| **Single-end**   | Each DNA fragment is sequenced only from one end.                           | RNA-seq with high depth, simple quantification. |
| **Paired-end**   | Both ends of a fragment are sequenced, providing more context and accuracy. | SNP discovery, structural variant detection. |
| **Long reads**   | Reads spanning thousands of bases, useful for repetitive or complex regions.| De novo genome assembly, isoform detection. |
| **HiFi reads**   | PacBio high-fidelity reads: long and highly accurate.                       | Reference-quality genomes, variant calling. |
| **Ultra-long reads** | Nanopore reads that can span >1 Mb.                                    | Telomere-to-telomere (T2T) assemblies.      |

---

### Quality Control Tools (e.g., FastQC, MultiQC, Trimmomatic)

Once raw sequencing data is generated, the first critical step is Quality Control (QC). QC assesses the quality of the raw reads and identifies potential issues that could affect downstream analysis. Poor quality data can lead to inaccurate alignments, incorrect variant calls, and misleading expression quantification.

**FastQC:**

FastQC is a widely used tool for performing quality control checks on raw sequencing data (FASTQ files). It generates a comprehensive report with various modules that highlight potential problems in your data, such as:

*   **Per base sequence quality:** Shows the quality scores across each position in the reads.
*   **Per sequence quality scores:** Distribution of average quality scores over all reads.
*   **Per base sequence content:** Checks for biases in nucleotide composition at each position.
*   **Per sequence GC content:** Distribution of GC content across all reads.
*   **Sequence Length Distribution:** Shows the distribution of read lengths.
*   **Overrepresented sequences:** Identifies sequences that appear unusually frequently, often indicating contamination or adapter sequences.
*   **Adapter content:** Detects the presence of common sequencing adapter sequences.

**MultiQC:**

MultiQC is a powerful tool that aggregates results from multiple bioinformatics analysis tools (including FastQC, but also aligners, variant callers, etc.) into a single, interactive HTML report. This makes it incredibly easy to get a quick overview of the quality and success of your entire workflow across many samples. Instead of opening dozens of individual FastQC reports, MultiQC consolidates them into one summary.

**Trimmomatic:**

Trimmomatic is a flexible and high-performance tool for trimming and cropping Illumina (and other) sequencing reads. It is used to remove low-quality bases, adapter sequences, and short reads that can negatively impact alignment and downstream analysis. Common trimming operations include:

*   **ILLUMINACLIP:** Removes adapter sequences.
*   **LEADING/TRAILING:** Removes low-quality bases from the start/end of reads.
*   **SLIDINGWINDOW:** Performs a sliding window quality trim.
*   **MINLEN:** Removes reads shorter than a specified length.

**Workflow Example: QC and Trimming**

1.  **Run FastQC on raw FASTQ files:**
    ```bash
    fastqc raw_reads_R1.fastq.gz raw_reads_R2.fastq.gz -o fastqc_reports
    ```
2.  **Run Trimmomatic to clean reads:**
    ```bash
    java -jar /path/to/trimmomatic.jar PE -phred33 \
        raw_reads_R1.fastq.gz raw_reads_R2.fastq.gz \
        trimmed_R1_paired.fastq.gz trimmed_R1_unpaired.fastq.gz \
        trimmed_R2_paired.fastq.gz trimmed_R2_unpaired.fastq.gz \
        ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    ```
    *Note: Replace `/path/to/trimmomatic.jar` and `adapters.fa` with actual paths.*

3.  **Run FastQC again on trimmed reads to assess improvement:**
    ```bash
    fastqc trimmed_R1_paired.fastq.gz trimmed_R2_paired.fastq.gz -o fastqc_trimmed_reports
    ```
4.  **Run MultiQC to summarize all FastQC reports:**
    ```bash
    multiqc fastqc_reports fastqc_trimmed_reports -o multiqc_summary
    ```

**Afternoon Session: Linux & Apocrita [Queen Mary’s HPC] Hands-on Challenge**

This session will be a practical challenge where you will apply your accumulated Linux command-line skills and get a taste of working in an HPC environment. We will use a simulated or actual HPC environment (like Queen Mary University of London\\'s Apocrita cluster, if access is provided and configured) to perform a basic bioinformatics task, emphasizing job submission and resource management.

**HPC Environment Overview (Apocrita Example):**

*   **Login:** Accessing the cluster via SSH.
*   **File Transfer:** Moving data to and from the cluster (e.g., using `scp` or `rsync`).
*   **Modules:** Loading necessary software (e.g., `module load FastQC`).
*   **Job Submission:** Writing and submitting a job script to the scheduler (e.g., Slurm).

**Example Slurm Job Script (`my_fastqc_job.sh`):**

```bash
#!/bin/bash
#SBATCH --job-name=fastqc_run      # Job name
#SBATCH --partition=long           # Partition name (e.g., short, long, gpu)
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks-per-node=1        # Number of tasks (cores) per node
#SBATCH --cpus-per-task=4          # Number of CPU cores per task
#SBATCH --mem=8G                   # Memory per node (e.g., 8GB)
#SBATCH --time=0-02:00:00          # Wall clock time limit (D-HH:MM:SS)
#SBATCH --output=fastqc_%j.out     # Standard output file
#SBATCH --error=fastqc_%j.err      # Standard error file

# Load necessary modules
module load fastqc/0.11.9 # Example version, check available modules

# Define input and output directories
INPUT_DIR="/path/to/your/raw_data"
OUTPUT_DIR="/path/to/your/fastqc_results"

mkdir -p $OUTPUT_DIR

# Run FastQC
fastqc ${INPUT_DIR}/sample_R1.fastq.gz ${INPUT_DIR}/sample_R2.fastq.gz -o $OUTPUT_DIR

echo "FastQC job completed!"
```

**Submitting the Job:**

```bash
ssh your_username@apocrita.qmul.ac.uk # Login to the HPC
sbatch my_fastqc_job.sh             # Submit the job
sq                                  # Check job status (Slurm queue)
```

**Challenge:**

Your challenge will be to take a provided raw sequencing dataset, perform quality control using FastQC and Trimmomatic, and then summarize the results using MultiQC, all within the HPC environment. This will test your ability to navigate the cluster, submit jobs, and manage your data effectively.

By the end of Week 1, you will have a strong foundation in the computational aspects of bioinformatics, including command-line proficiency, basic scripting, understanding of biological sequences and databases, and critical skills in data quality control and HPC usage. This knowledge will be essential as we move into more applied genomic analyses in Week 2.


## You have completed **Day 4**!

### Useful Database Tutorials

- [NCBI Databases Overview](https://www.ncbi.nlm.nih.gov/guide/all/)  
- [Entrez Direct User Guide](https://www.ncbi.nlm.nih.gov/books/NBK179288/)  
- [EBI Databases](https://www.ebi.ac.uk/services/)
- [Ensembl Genome Browser](https://www.ensembl.org/index.html)  
