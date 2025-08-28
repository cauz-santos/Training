## Week 1: Foundations of Bioinformatics and Genomic Data Handling

## Day 4: Experimental Design, Sequencing and Data QC 

Understanding how sequencing experiments are designed, executed, and validated is fundamental in bioinformatics.  
Poor experimental design or low-quality data can undermine even the most advanced analyses, leading to biased results or false discoveries.  

In this session, we will:  
- Explore the principles of **experimental design** and their impact on SNP and RNA-seq studies.  
- Compare different **sequencing technologies** and **read types**, discussing their strengths and limitations.  
- Learn about essential **quality control tools** (FastQC, Trimmomatic, MultiQC) to ensure reliable data.  
- Gain hands-on experience running QC workflows in a **High-Performance Computing (HPC)** environment.  

By the end of Day 4, you will understand not just **how sequencing works**, but also how to critically assess the **quality and reliability of sequencing data** — a crucial step for all downstream genomic analyses.

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

---
## Exercises: Sequencing Data Analysis  

This repository contains practical exercises to explore and analyze sequencing data, focusing on FASTQ files and UNIX command-line tools.


### Section 1: Exploring FASTQ Files with UNIX Commands  

In this section, you will learn to use basic UNIX commands to inspect and understand the FASTQ format, as well as perform preliminary quality analyses.

Navigate to your training data folder:

```bash
cd /media/uni/data2/Training/data/RADseq
ls -lh
```


**Exercise 1.1: Inspect a FASTQ file**  
Use the `less` command (or `zless` for compressed files) to inspect a FASTQ file. Observe the file structure and identify the initial header characters for each read.


```bash
zless DRR070477.fastq.gz
```

**Questions:**
*   What character does each read start with?
*   How many lines per read?

Remember:
FASTQ files are the standard format for storing raw sequencing data. Each read in a FASTQ file consists of four lines:

1.    **Sequence Identifier**: Starts with @ and includes metadata.  
2.    **Nucleotide Sequence**: The sequence of DNA bases (A, T, G, C).  
3.    **Separator Line**: Starts with + and may repeat the identifier.  
4.    **Quality Scores**: ASCII-encoded quality values corresponding to each base in the sequence.  
![image](https://github.com/user-attachments/assets/5ffbd3ff-adb8-4fb0-81bc-4c66609771bd)


**Exercise 1.2: Count number of reads**  
Determine the number of reads in a FASTQ file using the `zgrep` command. The command below counts lines that start with the specified header characters (`@SR`).

```bash
zgrep -c "^@SR" DRR070477.fastq.gz 
```


**Exercise 1.3: Check read length distribution**  
This `awk` command processes a compressed FASTQ file to analyze the distribution of read lengths. It first decompresses the file using `zcat`, then uses `awk` to extract the sequence lines (every 4th line, starting with the second). For each sequence line, it calculates the length of the sequence. The resulting lengths are then sorted and counted using `sort` and `uniq -c`, respectively, providing a summary of how many reads have each specific length.

```bash
zcat DRR070477.fastq.gz | awk '(NR%4==2){print length($0)}' | sort | uniq -c
```

**Exercise 1.4: Getting familiar with quality values** 

Quality values in FASTQ files are encoded using the Phred format. To convert an ASCII quality character to a numerical Phred value, you can use the ASCII table and subtract 33 (for Phred+33).

Use the ASCII table below to get the numerical value for the ASCII character `*`.

![image](https://github.com/user-attachments/assets/d86910fc-c7f7-4208-b8d4-9991678686f9)

It is 42. The convention is to subtract 33 (phred33 encoded), which makes 9. Is this a good quality? To find out the `p` value, calculate:

![image](https://github.com/user-attachments/assets/43dab0bc-eda2-4ec2-9ec0-dc5498ea4c4d)

Now that you know how to convert quality values, fill out the following table:

| Quality in fastq | Q in decimal | p |
|---|---|---|
| * | 9 | 0.1259 |
| I | 40 | |


____
### Section 2: Run FastQC Locally  

Even though we are on the cluster, this step is very fast and does not need to be submitted with Slurm.  
We just need to **load the module** and run FastQC directly.  

```bash
module load fastqc

mkdir -p fastqc_reports
fastqc DRR070477.fastq.gz -o fastqc_reports
```

This will create two output files inside fastqc_reports/:  
`DRR070477_fastqc.html` → open in a browser  
`DRR070477_fastqc.zip` → contains all raw QC data  

**Viewing the results on your Laptop**  
The cluster does not have a graphical interface, so you cannot open the .html report directly there.
To inspect the results:

On your laptop, create a folder for today’s results (e.g., inside your course folder):

```bash
mkdir -p ~/bioinformatics_training/day4
cd ~/bioinformatics_training/day4
```

**Copy the results from the cluster** to your local folder using `scp:

```bash
scp your_username@login02.lisc.univie.ac.at:/path/to/fastqc_reports/DRR070477_fastqc.html .
scp your_username@login02.lisc.univie.ac.at:/path/to/fastqc_reports/DRR070477_fastqc.zip .
```

Replace:  
`your_username` with your cluster username.  
`/path/to/fastqc_reports/` with the actual path where the results were created on the cluster.  
You can find the correct path by navigating to the folder in the cluster and typing:  

```bash
pwd
```

Open the .html file by double-clicking it or dragging it into your browser.

**Questions to Guide FastQC Interpretation:**

**Basic Statistics:**
*   How many sequences are there?
*   What is the sequence length?
*   What is the GC content?

**Per base sequence quality / Per sequence quality scores:**
*   Are there red ❌ or orange ⚠️ warnings?
*   Do you think this run produced good quality sequences?

**Per base sequence content / GC content:**
*   Do you observe any unusual base composition bias?
*   Should we worry about this in this particular case?

**Overrepresented sequences:**
*   Why does FastQC give a warning message here?
*   Could this be related to adapters or PCR duplicates?

___
### Section 3: Moving to the Cluster (Slurm Jobs with For Loops)

Now we will repeat the steps on the cluster using Slurm. Instead of running files one by one, we will use `for` loops to process them all.

**Step 3.1 – Run FastQC on all files**  
Open a new file called `fastqc_job.sh` with the text editor `vi`:

```bash
vi fastqc_job.sh
```

When the file opens, you are in command mode.
To start typing, press the `i` key (this puts you in insert mode).

Now copy and paste the following script into the file:

```bash
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH -o fastqc.out
#SBATCH -e fastqc.err

module load fastqc

mkdir -p fastqc_reports

for fq in *.fastq.gz
do
    echo "Running FastQC on $fq"
    fastqc "$fq" -o fastqc_reports
done
```

Save and exit (`ESC`, then type `:wq` and press `ENTER`).

Submit the job typing:

```bash
sbatch fastqc_job.sh
```

This will run FastQC on all `.fastq.gz` files in the current directory and save the output reports (.html and .zip) in the folder `fastqc_reports/.`



**Step 3.2 – View results on your laptop**  
Just like in Section 2, copy the results to your computer:

```bash
scp your_username@login02.lisc.univie.ac.at:/path/to/fastqc_reports/*.html ~/bioinformatics_training/day4/
```

Then open them in your browser to inspect the results.



**Step 3.3 – Trim Reads with Trimmomatic**  
Trimmomatic is a fast, multithreaded tool that removes adapters, trims poor-quality bases, and filters short reads. It works with both single-end and paired-end reads, and can handle compressed files (.gz).

**Common Trimmomatic Parameters:**

| Parameter | Description |
|---|---|
| `ILLUMINACLIP` | Cuts adapter and Illumina-specific sequences from the read. |
| `SLIDINGWINDOW` | Trims once the average quality within a sliding window falls below a threshold. |
| `MAXINFO` | Adaptive trimming balancing read length and error rate. |
| `LEADING` | Cuts bases from the start if below a threshold quality. |
| `TRAILING` | Cuts bases from the end if below a threshold quality. |
| `CROP` | Cuts the read to a fixed length by trimming from the end. |
| `HEADCROP` | Cuts a specified number of bases from the start. |
| `MINLEN` | Discards the read if shorter than the threshold length. |
| `AVGQUAL` | Drops the read if its average quality is below the threshold. |

Create a file `trim_job.sh`, with the text editor `vi`:

```bash
vi trim_job.sh
```

When the file opens, you are in command mode.
To start typing, press the `i` key (this puts you in insert mode).

Now copy and paste the following script into the file::

```bash
#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:30:00
#SBATCH -o trim.out
#SBATCH -e trim.err

module load trimmomatic

for fq in *.fastq.gz
do
    base=$(basename "$fq" .fastq.gz)
    echo "Trimming $fq → ${base}_trimmed.fastq.gz"
    java -jar $TRIMMOMATIC_JAR SE -phred33 \
      "$fq" "${base}_trimmed.fastq.gz" \
      ILLUMINACLIP:adapters.fa:2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
done
```

Submit the job:

```bash
sbatch trim_job.sh
```

**Step 2.3.3 – Run FastQC again on trimmed data**

Create `fastqc_trimmed_job.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=fastqc_trimmed
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH -o fastqc_trimmed.out
#SBATCH -e fastqc_trimmed.err

module load fastqc

mkdir -p fastqc_trimmed_reports

for fq in *_trimmed.fastq.gz
do
    echo "Running FastQC on $fq"
    fastqc "$fq" -o fastqc_trimmed_reports
done
```

Submit the job:

```bash
sbatch fastqc_trimmed_job.sh
```

Inspect the results:

Download the reports to your local machine (as in Section 2):

```bash
scp your_username@login02.lisc.univie.ac.at:/path/to/fastqc_trimmed_reports/*.html ~/bioinformatics_training/day4/
```

Open the `.html` files in your browser.

**Step 2.3.4 – Summarize with MultiQC**

Create `multiqc_job.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH -o multiqc.out
#SBATCH -e multiqc.err

module load multiqc/1.14

multiqc fastqc_reports fastqc_trimmed_reports -o multiqc_summary
```

Submit the job:

```bash
sbatch multiqc_job.sh
```

### 2.4: Reflection Questions

*   How many reads were removed after trimming?
*   Did the average base quality improve?
*   Were adapters successfully removed?
*   Why does read length distribution change after trimming?
*   Why is QC important before mapping or SNP calling?

**Challenge:**

Your challenge will be to take a provided raw sequencing dataset, perform quality control using FastQC and Trimmomatic, and then summarize the results using MultiQC, all within the HPC environment. This will test your ability to navigate the cluster, submit jobs, and manage your data effectively.

By the end of Week 1, you will have a strong foundation in the computational aspects of bioinformatics, including command-line proficiency, basic scripting, understanding of biological sequences and databases, and critical skills in data quality control and HPC usage. This knowledge will be essential as we move into more applied genomic analyses in Week 2.


## You have completed **Day 4**!

### Useful Database Tutorials

- [NCBI Databases Overview](https://www.ncbi.nlm.nih.gov/guide/all/)  
- [Entrez Direct User Guide](https://www.ncbi.nlm.nih.gov/books/NBK179288/)  
- [EBI Databases](https://www.ebi.ac.uk/services/)
- [Ensembl Genome Browser](https://www.ensembl.org/index.html)  
