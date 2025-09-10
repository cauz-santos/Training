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

**Table 1. Types of Sequencing and Reads**

| Sequencing Type        | Description                                                                 | Typical Applications                          | Read Types & Lengths                |
|------------------------|-----------------------------------------------------------------------------|-----------------------------------------------|--------------------------------------|
| **Sanger Sequencing**  | First-generation sequencing, low throughput but very accurate.              | Small gene fragments, validation, cloning.    | Single reads, ~500–1000 bp           |
| **Illumina (NGS)**     | Short-read sequencing with high accuracy and throughput.                    | WGS, RNA-seq, SNP genotyping, metagenomics.   | Paired-end (2 × 50–300 bp) or single-end |
| **Ion Torrent**        | Short-read sequencing, detects pH changes during DNA synthesis.             | Small genomes, targeted sequencing.           | Single or paired-end, up to ~400 bp  |
| **PacBio (HiFi/CLR)**  | Long-read sequencing; HiFi reads combine length with high accuracy.         | Genome assembly, isoform sequencing, SVs.     | Continuous Long Reads (up to 50 kb) or HiFi (~10–25 kb, >99% accurate) |
| **Oxford Nanopore (ONT)** | Portable, real-time long-read sequencing.                               | Structural variants, genome assembly, direct RNA sequencing. | Reads from kb to Mb (very long), lower raw accuracy but improving with new chemistries |
| **10x Genomics (Linked-Reads)** | Uses barcoding to connect short Illumina reads to long molecules. | Phasing, haplotyping, structural variation.   | Short reads (~150 bp) but linked across long DNA molecules (~50–100 kb) |


**Table 2. Types of Reads**

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

First create a folder for the file outputs of day 6 in your home directory:
   ```bash
   mkdir 04_qc_trimming
   ```
Now, enter the folder with the command `cd`

Then create some subfolders, for each specific analysis we will perform:
   ```bash
   mkdir fastq trimmed
   ```

### Section 1: Exploring FASTQ Files with UNIX Commands  

In this section, you will learn to use basic UNIX commands to inspect and understand the FASTQ format, as well as perform preliminary quality analyses.

Navigate to our training data folder:

```bash
cd /lisc/scratch/course/pgbiow/data/RADseq
ls -lh
```


**Exercise 1.1: Inspect a FASTQ file**  
Use the `less` command (or `zless` for compressed files) to inspect a FASTQ file. Observe the file structure and identify the initial header characters for each read.


```bash
zless EO_Ind1.fastq.gz
```

**❓ Questions:**
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
zgrep -c "^@DR" EO_Ind1.fastq.gz
 
```

Press `q` to exit

**Exercise 1.3: Check read length distribution**  
This `awk` command processes a compressed FASTQ file to analyze the distribution of read lengths. It first decompresses the file using `zcat`, then uses `awk` to extract the sequence lines (every 4th line, starting with the second). For each sequence line, it calculates the length of the sequence. The resulting lengths are then sorted and counted using `sort` and `uniq -c`, respectively, providing a summary of how many reads have each specific length.

```bash
zcat EO_Ind1.fastq.gz | awk '(NR%4==2){print length($0)}' | sort | uniq -c
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
We just need to **load the module**, change the output directory to your home directory  and run FastQC directly.  

```bash
module load fastqc


fastqc EO_Ind1.fastq.gz -o /path/to/your/home/directory/04_qc_trimming/fastq/
```

This will create two output files inside fastqc_reports/:  
`EO_Ind1_fastqc.html` → open in a browser  
`EO_Ind1_fastqc.zip` → contains all raw QC data  

**Viewing the results on your Laptop**  
The cluster does not have a graphical interface, so you cannot open the .html report directly there.
To inspect the results:

On your laptop, create a folder for today’s results (e.g., inside your course folder):

```bash
mkdir -p ~/bioinformatics_training/day4/fastqc
cd ~/bioinformatics_training/day4/fastqc
```

**Copy the results from the cluster** to your local folder using `scp:

```bash
scp your_username@login02.lisc.univie.ac.at:/path/to/fastqc_reports/EO_Ind1_fastqc.html .
scp your_username@login02.lisc.univie.ac.at:/path/to/fastqc_reports/EO_Ind1_fastqc.zip .
```

Replace:  
`your_username` with your cluster username.  
`/path/to/fastqc_reports/` with the actual path where the results were created on the cluster.  
You can find the correct path by navigating to the folder in the cluster and typing:  

```bash
pwd
```

Open the .html file by double-clicking it or dragging it into your browser.

**❓ Questions to Guide FastQC Interpretation:**

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
### Section 3: Moving to Slurm Jobs with For Loops

Now we will repeat the steps on the cluster using Slurm. Instead of running files one by one, we will use `for` loops to process them all.

Please move back to the folder `04_qc_trimming` in your home directory:

```bash
cd path/to/your/home/directory/04_qc_trimming/
```

**Step 3.0 – Instructions on Creating and Editing Shell Scripts with `vi`**

Throughout this training, we are creating and editing shell scripts. We will use the `vi` editor as in the previous days, a powerful text editor available on most UNIX-like systems. Here's a quick guide:

1.  **Open or Create a File:**
    To create a new file or open an existing one, type:
    ```bash
    vi <filename>
    ```
    For example, to create `index_genome.sh`:
    ```bash
    vi fastqc_job.sh
    ```

2.  **Insert Mode:**
    When `vi` opens, you are in **command mode**. To start typing, you need to enter **insert mode**. Press the `i` key.
    You should see `-- INSERT --` at the bottom of your terminal.

3.  **Type your content:**
    Now you can type or paste the script content.

4.  **Exit Insert Mode:**
    Once you are done typing, press the `Esc` key to return to **command mode**.

5.  **Save and Quit:**
    In command mode, type `:wq` (write and quit) and press `Enter`.
    *   `:w` saves the file.
    *   `:q` quits `vi`.
    *   `:wq` saves and quits.
    *   `:q!` quits without saving (use with caution!).



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

# Define input and output directories
INPUT_DIR="/lisc/scratch/course/pgbiow/data/RADseq"
OUTPUT_DIR="/lisc/scratch/course/pgbiow/04_qc_trimming/fastq"

# Create output directory if it doesn’t exist
mkdir -p "$OUTPUT_DIR"

# Loop through all FASTQ files in the input directory
for fq in "$INPUT_DIR"/*.fastq.gz
do
    echo "Running FastQC on $fq"
    fastqc "$fq" -o "$OUTPUT_DIR"
done
```

Save and exit (`ESC`, then type `:wq` and press `ENTER`).

Submit the job typing:

```bash
sbatch fastqc_job.sh
```

You can check how the analysis is progressing in the `fastqc.err` file:
```bash
sbatch fastqc.err
```

This will run FastQC on all `.fastq.gz` files in the current directory and save the output reports (.html and .zip) in the folder `fastqc`



**Step 3.2 – View results on your laptop**  
Just like in Section 2, copy the results to your computer:

```bash
scp your_username@login02.lisc.univie.ac.at:/path/to/fastqc_reports/*.html ~/bioinformatics_training/day4/fastqc
```

Then open them in your browser to inspect the results.

**Step 3.3 – Summarize with MultiQC**  
After running FastQC on both **raw reads** and **trimmed reads**, you now have multiple `.html` reports (one per file).  
Opening them one by one is time-consuming and makes it difficult to compare.  

**MultiQC** solves this problem by **aggregating all FastQC reports** into a single, easy-to-navigate HTML summary.  
This allows you to directly compare quality improvements before and after trimming.

Create a new file called `multiqc_job_raw.sh`:

```bash
vi multiqc_job_raw.sh
```

Press `i` to enter insert mode and paste the script:

```bash
#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH -o multiqc.out
#SBATCH -e multiqc.err

set -euo pipefail

# conda for non-interactive jobs
source /lisc/app/conda/miniforge3/etc/profile.d/conda.sh

ENV_NAME="multiqc-1.30"
RAW_DIR="/lisc/scratch/course/pgbiow/04_qc_trimming/fastq"
OUT_DIR="/lisc/scratch/course/pgbiow/04_qc_trimming/multiqc_raw"

mkdir -p "$OUT_DIR"

# Optional sanity check
conda run --no-capture-output -n "$ENV_NAME" multiqc --version

# Run on the raw FastQC outputs in RAW_DIR
conda run --no-capture-output -n "$ENV_NAME" multiqc "$RAW_DIR" \
  -o "$OUT_DIR" \
  -n "multiqc_fastqc_raw_$(date +%F).html" \
  --force
```

Save and exit (`ESC`, then `:wq`).

Submit the job:

```bash
sbatch multiqc_job_raw.sh
```

**Step 3.4 – Copy the MultiQC report to your local computer**  
On your laptop, create a folder for results:

```bash
mkdir -p ~/bioinformatics_training/day4/multiqc
cd ~/bioinformatics_training/day4/multiqc
```

Then copy the HTML file from the cluster:

```bash
scp your_username@login02.lisc.univie.ac.at:/path/to/multiqc_summary/*.html .
```

Replace `/path/to/multiqc_summary/` with the actual path on the cluster for the `multiqc_raw` folder (check with pwd).

Finally, open the report in your web browser by double-clicking it.


**Step 3.5 – Trim Reads with Trimmomatic**  
Trimmomatic is a fast, multithreaded tool that removes adapters, trims poor-quality bases, and filters short reads. It works with both single-end and paired-end reads, and can handle compressed files (.gz).

**Table 3. Common Trimmomatic Parameters**

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
#SBATCH --job-name=trimmomatic_strict
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:30:00
#SBATCH -o trim.out
#SBATCH -e trim.err

module load trimmomatic

INPUT_DIR="/lisc/scratch/course/pgbiow/data/RADseq"
OUTPUT_DIR="/path/to/your/home/directory/04_qc_trimming/trimmed"

TRIMMOMATIC_JAR="/lisc/app/trimmomatic/0.39/trimmomatic-0.39.jar"
ADAPTERS="/lisc/app/trimmomatic/0.39/adapters/TruSeq3-SE.fa"

mkdir -p "$OUTPUT_DIR"

for fq in "$INPUT_DIR"/*.fastq.gz
do
    base=$(basename "$fq" .fastq.gz)
    out="${OUTPUT_DIR}/${base}_trimmed.fastq.gz"
    log="${OUTPUT_DIR}/${base}.trimlog.txt"

    echo "Trimming $fq → $out"

    java -jar "$TRIMMOMATIC_JAR" SE \
      -threads "$SLURM_CPUS_PER_TASK" -phred33 -trimlog "$log" \
      "$fq" "$out" \
      ILLUMINACLIP:${ADAPTERS}:1:30:15:8 \
      LEADING:20 TRAILING:20 \
      SLIDINGWINDOW:4:25 \
      AVGQUAL:25 \
      MINLEN:50
done
```

Submit the job:

```bash
sbatch trim_job.sh
```

**Step 3.6 – Checking Trimmomatic Results**  
After running the trimming job, Trimmomatic reports how many reads were processed, how many survived, and how many were dropped. These statistics are written to the log file `trim.err`.


**Inspect the trimming log:**  
Check the log file to see the output for each sample:

```bash
less trim.err
```

Look for lines like:
```bash
Input Reads: 1871406 Surviving: 1682221 (89.89%) Dropped: 189185 (10.11%)
```

**Extract only the summary lines:**  
To quickly show all trimming results:
```bash
grep "Input Reads" trim.err
```

**Create a clean table summary:**  
We can reformat the results into a simple table using awk and save it to a file:
```bash
grep "Input Reads" trim.err | \
awk '{print "Sample " NR ": Input="$3", Surviving="$5", Dropped="$7}' \
> trimming_summary.txt
```

**View the summary file:**  
Finally, display the table:
```bash
cat trimming_summary.txt
```

Example output:
```bash
Sample 1: Input=1871406, Surviving=1682221, Dropped=189185
Sample 2: Input=1923400, Surviving=1749287, Dropped=174113
```


**Step 3.7 – Run FastQC again on trimmed data**  
Create `fastqc_trimmed_job.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=fastqc_trimmed
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH -o fastqc_trimmed.out
#SBATCH -e fastqc_trimmed.err

module load fastqc

# Input and output directories
INPUT_DIR="/lisc/scratch/course/pgbiow/04_qc_trimming/trimmed"
OUTPUT_DIR="/lisc/scratch/course/pgbiow/04_qc_trimming/fastqc_trimmed_reports"

# Create output directory if it doesn’t exist
mkdir -p "$OUTPUT_DIR"

# Loop through all trimmed FASTQ files
for fq in "$INPUT_DIR"/*_trimmed.fastq.gz
do
    echo "Running FastQC on $fq"
    fastqc "$fq" -o "$OUTPUT_DIR"
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

**Step 3.8 – Summarize with MultiQC**  
Create a new file called `multiqc_job_trimmed.sh`:

```bash
vi multiqc_job_trimmed.sh
```

Press `i` to enter insert mode and paste the script:

```bash
#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH -o multiqc.out
#SBATCH -e multiqc.err

set -euo pipefail

# conda for non-interactive jobs
source /lisc/app/conda/miniforge3/etc/profile.d/conda.sh

ENV_NAME="multiqc-1.30"
RAW_DIR="/lisc/scratch/course/pgbiow/04_qc_trimming/trimmed"
OUT_DIR="/lisc/scratch/course/pgbiow/04_qc_trimming/multiqc_trimmed"

mkdir -p "$OUT_DIR"

# Optional sanity check
conda run --no-capture-output -n "$ENV_NAME" multiqc --version

# Run on the raw FastQC outputs in RAW_DIR
conda run --no-capture-output -n "$ENV_NAME" multiqc "$RAW_DIR" \
  -o "$OUT_DIR" \
  -n "multiqc_fastqc_raw_$(date +%F).html" \
  --force
```

Save and exit (`ESC`, then `:wq`).

Submit the job:

```bash
sbatch multiqc_job_trimmed.sh
```

**Step 3.9 – Copy the MultiQC report to your local computer**  
On your laptop, copy the HTML file from the cluster:

```bash
scp your_username@login02.lisc.univie.ac.at:/path/to/multiqc_summary/*.html .
```

Replace `/path/to/multiqc_summary/` with the actual path on the cluster for the `multiqc_trimmed` folder (check with pwd).

Finally, open the report in your web browser by double-clicking it.


### ❓ Questions:

*   How many reads were removed after trimming? Was trimming too aggressive (discarding too many reads)? How would you adjust parameters (e.g., SLIDINGWINDOW or MINLEN) if needed?
*   Did the average base quality changed?
*   Were adapters successfully removed?
*   Why does read length distribution change after trimming?
*   Why is QC important before mapping or SNP calling?


## Final Tips for Sequencing Data QC

- **Always keep raw data untouched**  
  Store original FASTQ files in a separate folder. Perform QC and trimming on *copies*, so you can always go back to the raw data.

- **Check sequencing depth early**  
  Knowing how many reads you have (and their length) helps decide if coverage is sufficient for your project (e.g., SNP discovery, RNA-seq, genome assembly).

- **Interpret FastQC, don’t just run it**  
  The red ❌ or orange ⚠️ are only flags. Ask yourself: *Does this affect my analysis?* (e.g., adapter contamination matters for SNP calling, but GC bias may not).

- **Trim carefully**  
  Over-trimming can discard too many reads, while under-trimming leaves poor-quality data. Adjust parameters like `SLIDINGWINDOW` and `MINLEN` to balance quality and read retention.

- **Always re-check after trimming**  
  Run FastQC again (or MultiQC for summaries) to confirm quality improvements and ensure adapters/low-quality bases were removed.

- **Use MultiQC for comparisons**  
  Instead of opening dozens of reports, use MultiQC to get a single overview. This is especially important in projects with many samples.

- **Document everything**  
  Keep a log of the commands, parameters, and versions of software used. This ensures reproducibility and helps troubleshoot issues later.

- **QC before analysis**  
  Poor QC can lead to false SNPs, biased expression levels, or failed assemblies. High-quality input is the foundation for reliable downstream analyses.

---
💡 **Remember:** Good experimental design and rigorous QC save time, money, and frustration in the long run. Treat this step as *essential*, not optional.


## You have completed **Day 4**!

###  Useful Tutorials on Sequencing & Data Quality

- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) – Official user guide with explanations of each module and how to interpret plots.  
- [MultiQC Documentation](https://multiqc.info/docs/) – Learn how to aggregate QC results from multiple samples and tools.  
- [Trimmomatic Manual](http://www.usadellab.org/cms/?page=trimmomatic) – Detailed explanation of trimming steps and parameters.
- [Galaxy Training: Quality Control (FastQC/Cutadapt/MultiQC)](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html) — hands-on QC workflow.  
- [HBC Training: QC of Raw NGS Data with FastQC](https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/lessons/QC_raw_data.html) — interpreting FastQC modules & flags.  
- [Wrangling Genomics: Trimming & Filtering](https://bookdown.org/ggiaever/wrangling-genomics/trimming-and-filtering.html) — didactic walkthrough of Trimmomatic options with examples.
 

