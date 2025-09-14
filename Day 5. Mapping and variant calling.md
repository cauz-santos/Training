## Week 1: Foundations of Bioinformatics and Genomic Data Handling

## Day 5: Genome mapping and variant calling 

Welcome to Day 5 of our bioinformatics training! Today, we'll delve into the core of genomic analysis: mapping sequencing reads to a reference genome and identifying genetic variations. This is a hands-on session designed for a high-performance computing (HPC) cluster environment using the Slurm workload manager.

### Learning Objectives

By the end of this session, you will be able to:

*   **Understand the principles** of read mapping and its importance.
*   **Prepare a reference genome** for alignment.
*   **Write and submit Slurm job scripts** to automate analysis on a cluster.
*   **Use `for` loops** to efficiently process multiple samples.
*   **Align sequencing reads** to a reference genome using BWA.
*   **Process and inspect alignment files** (SAM/BAM) using Samtools.
*   **Assess mapping quality** by checking mapping rates and coverage.
*   **Perform variant calling** using BCFtools to identify SNPs and indels.
*   **Analyze and filter VCF files** to obtain high-confidence variant calls.

### Session Overview

We will be working with a set of up to 13 single-end FASTQ files and a reference genome. The workflow is as follows:

1.  **Setup:** Prepare the environment and data.
2.  **Reference Genome Indexing:** Create an index for the reference genome to enable fast alignment.
3.  **Read Mapping:** Align the FASTQ reads to the indexed reference genome.
4.  **Alignment Post-processing:** Convert SAM to BAM, sort, and index the alignment files.
5.  **Quality Control of Alignments:** Check mapping rates and coverage.
6.  **Variant Calling:** Identify genetic variants from the aligned reads.
7.  **Variant Filtering and Annotation:** Filter raw variants to improve accuracy and interpret the results.

### The Data

*   **Reference Genome of *Elaeis guineensis*:** Elaeis_guineensis_genomic.fna` (Size: 1.8 GB)
*   **Input Reads:** 13 single-end FASTQ files trimmed in the previous day

**Table 1. File Format Flow**
| File Type | Extension | Description | Produced By |
|-----------|-----------|-------------|-------------|
| **Raw Reads** | `.fastq.gz` | Raw sequences with qualities | Sequencer |
| **Alignments (text)** | `.sam` | Large, human-readable alignments | BWA |
| **Alignments (binary)** | `.bam` | Compressed, efficient alignments | Samtools |
| **Sorted + Index** | `.sorted.bam` + `.bai` | Required for downstream analysis | Samtools |
| **Variants** | `.vcf.gz` | SNPs and indels compared to reference | BCFtools |

---

### Section 1: Setting up the Environment and Data

Before we start, let's make sure our environment is set up correctly on the HPC cluster.

1.  **Log in to the cluster** and navigate to your working directory.

2.  **Create a directory for Day 5:**

First create a folder for the file outputs of day 6 in your home directory:
   ```bash
   mkdir 05_mapping_varriant_calling
   ```
Now, enter the folder with the command `cd`

Then create some subfolders, for each specific analysis we will perform:
   ```bash
   mkdir reference bwa_mapping bcftools_variants
   ```


3.  **Load necessary modules:**

    We will need `bwa`, `samtools`, and `bcftools` for this session. The specific module names might vary on your cluster.

    ```bash
    module load bwa
    module load samtools
    module load bcftools
    ```

    **❓ Question:** Why do we use `module load`? What happens if we don't load the modules?

4.  **Data:** Ensure your FASTQ files and the reference genome are in your working directory or accessible via a path.


### Section 2: Reference Genome Indexing

**Why do we index the reference genome?**

Aligning reads to a large reference genome (1.8 GB in our case) is computationally intensive. Searching for the best match for millions of short reads against the entire genome sequence would be incredibly slow. Indexing creates a set of auxiliary files that act like a 'phonebook' for the genome, allowing the aligner (BWA) to quickly find the location of a specific sequence. This is a one-time step for each reference genome.

We will create a Slurm job script to perform the indexing.

**Creating and Editing Shell Scripts with `vi`**

Throughout this training, we are creating and editing shell scripts. We will use the `vi` editor as in the previous days, a powerful text editor available on most UNIX-like systems. Here's a quick guide:

1.  **Open or Create a File:**
    To create a new file or open an existing one, type:
    ```bash
    vi <filename>
    ```
    For example, to create `index_genome.sh`:
    ```bash
    vi index_genome.sh
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

Now, move to your home folder and enter:

```bash
cd 05_mapping_varriant_calling
```

let's create the `index_genome.sh` file:
```bash
#!/bin/bash
#SBATCH --job-name=bwa_index
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH -o bwa_index.out
#SBATCH -e bwa_index.err

# Load the BWA module
module load bwa

# Define input and output directories
INPUT_DIR="/lisc/scratch/course/pgbiow/data/genomes"
OUTPUT_DIR="/path/to/your/home/directory/05_mapping_varriant_calling/reference"

# Define the reference genome file name
GENOME="Elaeis_guineensis_genomic.fna"

# Create output directory if it doesn’t exist
mkdir -p "$OUTPUT_DIR"

# Copy the reference genome to the output directory
echo "Copying reference genome to $OUTPUT_DIR"
cp "$INPUT_DIR/$GENOME" "$OUTPUT_DIR/"

# Change to the output directory
cd "$OUTPUT_DIR" || exit 1

# Run BWA index on the copied genome
echo "Indexing the reference genome: $GENOME"
bwa index "$GENOME"

echo "Indexing complete."
```

**Submit the job:**

```bash
sbatch index_genome.sh
```

**❓ Question:** What do the `#SBATCH` directives in the script mean? Why are they important for running jobs on a cluster?

Once the job is complete, you will see several new files in your directory with extensions like `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`. These are the BWA index files.

___

### Section 3: Read Mapping and Post-processing


Now that we have our indexed genome, we can align our FASTQ reads.  
This step takes the **raw reads** and transforms them into **clean, analysis-ready alignment files**.  

The main steps are:  

1. **Read Mapping (BWA-MEM)** → aligns each read to the reference genome.  
2. **SAM → BAM conversion (Samtools view)** → compresses the large SAM text file into a smaller binary BAM file.  
3. **Sorting (Samtools sort)** → orders the reads by their position in the genome (required for many tools).  
4. **Duplicate Removal (Samtools markdup)** → removes PCR duplicates that can bias variant calling.  
   - We use the flag `-r` to **remove** duplicates instead of just marking them.  
5. **Indexing (Samtools index)** → creates a `.bai` index file so tools can quickly access specific regions of the BAM.  
6. **Cleanup** → intermediate files (unsorted BAMs) are deleted to save space.

**⚠️ Important note**: Duplicate removal is useful for SNP/variant calling.
But in other applications (e.g., RNA-seq or ChIP-seq), duplicates may carry biological information and should not be removed.


**Create a file named `map_reads.sh` using `vi`:**

```bash
#!/bin/bash
#SBATCH --job-name=bwa_picard_pipeline
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --partition=basic
#SBATCH --time=06:00:00
#SBATCH -o bwa_picard_pipeline.out
#SBATCH -e bwa_picard_pipeline.err

set -euo pipefail

# Load required modules
module purge
module load bwa
module load samtools
module load picard/3.4.0 || true         # may only provide the JAR
module load java/17 || module load openjdk/17 || true  # Picard 3.x wants Java 17

# Resolve Picard command (wrapper if present, else JAR)
if command -v picard &>/dev/null; then
  PICARD_CMD="picard"
else
  PICARD_JAR="/lisc/app/picard/3.4.0/picard.jar"
  if [ ! -f "$PICARD_JAR" ]; then
    echo "ERROR: Picard jar not found at $PICARD_JAR and no 'picard' wrapper in PATH." >&2
    exit 1
  fi
  # Heap ~50G; adjust if needed
  PICARD_CMD="java -Xmx50g -jar $PICARD_JAR"
fi

# Define directories
TRIMMED_DIR="/lisc/scratch/course/pgbiow/04_qc_trimming/trimmed"
REF_DIR="/lisc/scratch/course/pgbiow/05_mapping_varriant_calling/reference"
OUT_DIR="/lisc/scratch/course/pgbiow/05_mapping_varriant_calling/bwa_mapping"

# Reference genome
GENOME="$REF_DIR/Elaeis_guineensis_genomic.fna"

# Create output directory if it doesn’t exist
mkdir -p "$OUT_DIR"

# Ensure reference is indexed (safe to re-run)
[ -f "${GENOME}.bwt" ] || bwa index "$GENOME"
[ -f "${GENOME}.fai" ] || samtools faidx "$GENOME"

# Move to trimmed FASTQ directory
cd "$TRIMMED_DIR" || exit 1

# Loop through all FASTQ files
for fq in *.fastq.gz; do
  base=$(basename "$fq" .fastq.gz)
  echo "Processing sample: $base"

  # Step 1: Align reads using BWA-MEM (single-end)
  bwa mem -t "$SLURM_CPUS_PER_TASK" "$GENOME" "$fq" > "$OUT_DIR/${base}.sam"

  # Step 2: Convert SAM to BAM
  samtools view -@ "$SLURM_CPUS_PER_TASK" -b "$OUT_DIR/${base}.sam" > "$OUT_DIR/${base}.bam"

  # Step 3: Sort BAM
  samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$OUT_DIR/${base}.sorted.bam" "$OUT_DIR/${base}.bam"

  # Step 4: Mark (or remove) duplicates with Picard
  # Tip: for many variant-calling workflows you MARK duplicates (REMOVE_DUPLICATES=false).
  $PICARD_CMD MarkDuplicates \
    I="$OUT_DIR/${base}.sorted.bam" \
    O="$OUT_DIR/${base}.sorted.dedup.bam" \
    M="$OUT_DIR/${base}.dup_metrics.txt" \
    REMOVE_DUPLICATES=false \
    ASSUME_SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR="${TMPDIR:-$OUT_DIR}"

  # Step 5: Index final BAM
  samtools index "$OUT_DIR/${base}.sorted.dedup.bam"

  echo "Finished processing: $base"
done

echo "All samples processed."
```

**Submit the job:**

```bash
sbatch map_reads.sh
```

**📍 Checkpoint:** After completion, you should see these files in the folder:
`sample.sorted.bam`
`sample.sorted.bam.bai`


**Quick Inspection of the SAM File**    
Our pipeline now produces both the human-readable SAM file and the processed BAM file.
The SAM (Sequence Alignment/Map) format is a plain-text, tab-delimited file that describes how each read aligns to the reference genome.

It has two main parts:
- **Header section** (lines starting with @) with reference and alignment metadata.
- **Alignment section** (one line per read) describing the mapping result.

To inspect one of the SAM files, we can simply use head:

```bash
#### Look at the first 30 lines of a SAM file
head -n 30 EO_Ind10_trimmed.sam
```

You should see a few header lines, followed by the first alignments.
Because SAM files are very large, we normally don’t keep them long-term — they are mainly for teaching, demonstration, or debugging.

**Example line of a mapped read:**  
```bash
SRR123456.1    99    NC_012345.1    10468    60    76M    =    10550    158    ACTG...    IIII...
```

Here:
- `QNAME` = SRR123456.1 → read name
- `FLAG` = 99 → paired, properly mapped (bitwise decoded)
- `RNAME` = NC_012345.1 → reference sequence
- `POS` = 10468 → starting coordinate
- `MAPQ` = 60 → high mapping quality
- `CIGAR` = 76M → 76 bases matched
- `SEQ` and `QUAL` → the actual read and base qualities

**Table 2. The SAM file format fields**

| Col | Field   | Description |
|-----|---------|-------------|
| 1   | **QNAME** | Query (read) name |
| 2   | **FLAG**  | Bitwise flag encoding read properties (mapped, paired, etc.) |
| 3   | **RNAME** | Reference sequence name (chromosome/contig) |
| 4   | **POS**   | 1-based leftmost position of the read on the reference |
| 5   | **MAPQ**  | Mapping quality (Phred-scaled) |
| 6   | **CIGAR** | Encodes alignment (matches, insertions, deletions) |
| 7   | **MRNM**  | Mate reference sequence name (`=` if same as RNAME) |
| 8   | **MPOS**  | Mate position |
| 9   | **ISIZE** | Inferred insert size |
| 10  | **SEQ**   | Read sequence |
| 11  | **QUAL**  | Quality string (ASCII-encoded Phred scores) |
| 12  | **OPT**   | Optional fields (e.g., NM:i:1 = edit distance) |


**❓ Questions:**

*   Why do we pipe the output of `bwa mem` directly to `samtools view`? What is the advantage of this approach?
*   What is the difference between a SAM and a BAM file? Why do we prefer BAM?
*   Why is sorting the BAM file a necessary step?



### Section 4: Alignment Quality Control

After mapping, it's crucial to assess the quality of our alignments. We'll check two key metrics: **mapping rate** and **coverage**.

**4.1: Checking Mapping Rate**  

The mapping rate tells us what percentage of our reads successfully aligned to the reference genome. A low mapping rate could indicate problems with the sequencing data, the reference genome, or contamination.

We can use `samtools flagstat` to get a summary of mapping statistics.

**Run `samtools flagstat` for each sample:**  

You can do this in a `for` loop using the following script :
create the script file:
```bash
vi check_mapping_rate.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=mapping_stats
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH -o mapping_stats.out
#SBATCH -e mapping_stats.err

module load samtools

# Define directory with BAMs
BAM_DIR="/lisc/scratch/course/pgbiow/05_mapping_varriant_calling/bwa_mapping"
OUT_FILE="$BAM_DIR/mapping_statistics.txt"

# Go to BAM directory
cd "$BAM_DIR" || exit 1

# Remove old stats file if it exists
rm -f "$OUT_FILE"

# Loop over all deduplicated BAMs
for bam in *.sorted.dedup.bam; do
    echo "Statistics for $bam:" >> "$OUT_FILE"
    samtools flagstat "$bam" >> "$OUT_FILE"
    echo "--------------------" >> "$OUT_FILE"
done

echo "All statistics saved to $OUT_FILE"
```

Submit the job:
```bash
sbatch check_mapping_rate.sh
```

When it finishes, check the output file:
```bash
cat /lisc/scratch/course/pgbiow/05_mapping_varriant_calling/bwa_mapping/mapping_statistics.txt
```

**Look for the line that says `... + ... mapped (...%)`.** This is your overall mapping rate.

**❓ Question:** What would you consider a good mapping rate? What might be the cause of a low mapping rate?


**4.2: Checking Coverage**  

Coverage (or read depth) is the number of reads that align to, or "cover," a specific position in the genome. Adequate coverage is essential for confident variant calling. Low coverage can lead to false negatives (missing real variants), while extremely high coverage might indicate PCR duplicates or repetitive regions.

We can use `samtools depth` or `samtools coverage` to assess coverage.

**Calculate average coverage across the genome:**

create the script file:
```bash
vi check_coverage.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=coverage_stats
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH -o coverage_stats.out
#SBATCH -e coverage_stats.err

module load samtools

# Directories
BAM_DIR="/lisc/scratch/course/pgbiow/05_mapping_varriant_calling/bwa_mapping"
OUT_FILE="$BAM_DIR/coverage_statistics.txt"

cd "$BAM_DIR" || exit 1
rm -f "$OUT_FILE"

# Loop through all deduplicated BAMs
for bam in *.sorted.bam; do
    echo "Calculating coverage for $bam ..."
    
    # Method 1: Average coverage using samtools depth + awk
    avg_cov=$(samtools depth "$bam" | awk '{sum+=$3} END { if (NR>0) print sum/NR; else print 0 }')
    
    # Method 2 (faster alternative): samtools coverage
    # avg_cov=$(samtools coverage "$bam" | awk 'NR==2 {print $7}')

    echo -e "$bam\tAverage coverage = $avg_cov" >> "$OUT_FILE"
done

echo "All coverage results saved to $OUT_FILE"
```
Submit the job:
```bash
sbatch check_coverage.sh
```

When it finishes, check the output file:
```bash
cat /lisc/scratch/course/pgbiow/05_mapping_varriant_calling/bwa_mapping/coverage_statistics.txt
```

**A more comprehensive view with `samtools coverage`:**

This command provides a summary of coverage statistics, including the percentage of the genome covered at different depths.

```bash
# Run on a single BAM file
samtools coverage EO_Ind10_trimmed.sorted.dedup.bam
```

**❓ Question:** What is the difference between breadth of coverage and depth of coverage? Why are both important?

___

### Optional Exercise: Inspecting Alignments with Qualimap  

While `samtools flagstat` and `samtools coverage` give good summary statistics, sometimes it’s useful to generate more detailed reports with graphs and summaries. For this, we can use **Qualimap**.  

Qualimap provides detailed statistics about the mapping quality, coverage distribution, GC bias, and more.  
This is optional, but a great way to visually confirm the quality of your alignments.  

**Step 1 – Create a Slurm job for Qualimap**  

Open a new file:

```bash
vi qualimap_job.sh
```

Press `i` to enter insert mode, then paste:

```bash
#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH -o qualimap.out
#SBATCH -e qualimap.err

module load qualimap

# Path to BAM directory
BAM_DIR="/lisc/scratch/course/pgbiow/05_mapping_varriant_calling/bwa_mapping"

# Example: run on one BAM file
BAM="$BAM_DIR/EO_Ind10_trimmed.sorted.bam"

# Run Qualimap bamqc
qualimap bamqc -bam "$BAM" -outdir "$BAM_DIR/qualimap_report_EO_Ind10" -nt 4
```

Save and exit (`ESC`, then `:wq`).

**Step 2 – Submit the job**  
```bash
sbatch qualimap_job.sh
```

This will create a folder qualimap_report_sample1/ containing:
`qualimapReport.html` → open in your browser (after copying it to your laptop with scp)


**Step 3 – Copy the report to your local computer**
On your laptop:

On your laptop, create a folder for today’s results (e.g., inside your course folder):

```bash
mkdir -p ~/bioinformatics_training/day5/qualimap
cd ~/bioinformatics_training/day4/qualimap
```

```bash
scp your_username@login02.lisc.univie.ac.at:/path/to/qualimap_report_sample1/qualimapReport.html ~/bioinformatics_training/day5/
```

Then open the HTML report in your browser.

___

### Section 5: Variant Calling

Now that we have high-quality alignments, we can proceed to variant calling. This process identifies positions where the aligned reads consistently differ from the reference genome.

We will use `bcftools mpileup` and `bcftools call`.

*   `bcftools mpileup`: Summarizes the base calls of aligned reads at each position in the genome.
*   `bcftools call`: Calls SNPs and indels from the `mpileup` output.

**Create a Slurm script named `variant_calling.sh` using `vi`:**

```bash
#!/bin/bash
#SBATCH --job-name=variant_calling
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH -o variant_calling.out
#SBATCH -e variant_calling.err

module load bcftools

# Reference genome
GENOME="/lisc/scratch/course/pgbiow/05_mapping_varriant_calling/reference/Elaeis_guineensis_genomic.fna"

# Directories
BAM_DIR="/lisc/scratch/course/pgbiow/05_mapping_varriant_calling/bwa_mapping"
OUT_DIR="/lisc/scratch/course/pgbiow/05_mapping_varriant_calling/bcftools_variants"

# Create output directory if it doesn’t exist
mkdir -p "$OUT_DIR"

# Move to BAM directory
cd "$BAM_DIR" || exit 1

# Create list of BAM files
ls *.sorted.dedup.bam > bamlist.txt

# Joint variant calling
bcftools mpileup -Ou -f "$GENOME" -b bamlist.txt | \
bcftools call -mv -Oz -o "$OUT_DIR/joint_variants.vcf.gz"

# Index the VCF
bcftools index "$OUT_DIR/joint_variants.vcf.gz"

echo "Variant calling complete. Results saved to:"
echo "$OUT_DIR/joint_variants.vcf.gz"
```

**Submit the job:**

```bash
sbatch variant_calling.sh
```

**❓ Question:** What do the `-mv` options in `bcftools call` stand for? Why might you choose to call only multiallelic sites or only biallelic sites?



### Section 6: VCF File Inspection and Statistics

A VCF (Variant Call Format) file contains information about the identified variants. It's important to understand its structure and to check some basic statistics.

**6.1: Inspecting a VCF file**

You can view the contents of a gzipped VCF file using `zless` or `bcftools view`.

```bash
# View the header and first few variant lines
bcftools view bcftools_variants/joint_variants.vcf.gz | head -n 70
```

**The VCF file has two main parts:**

*   **Header:** Lines starting with `##` describe the file format, the reference genome, and the annotations used in the file (INFO, FORMAT fields).
*   **Variant entries:** Each line represents a variant, with columns for chromosome, position, reference allele, alternate allele, quality, and sample-specific information.

**❓ Question:** What is the difference between the `INFO` and `FORMAT` fields in a VCF file?

**6.2: VCF Statistics**

`bcftools stats` provides a wealth of information about your VCF file, including the number of SNPs and indels, transition/transversion ratio, and more.

```bash
bcftools stats bcftools_variants/joint_variants.vcf.gz > bcftools_variants/joint_variants.vcf.stats
```

You can then view the generated statistics file:

```bash
less bcftools_variants/joint_variants.vcf.stats
```

**❓ Question:** What is the transition/transversion (Ts/Tv) ratio, and why is it a useful quality metric for SNP calls?

**6.3: Checking Missing Data (per sample)**

Missing genotypes can bias downstream population analyses. Let’s quantify **percent missing per sample** on the *joint VCF*.

We’ll use **VCFtools**

```bash
#First load the module
module load vcftools

# Compute per-individual missingness
vcftools --gzvcf bcftools_variants/joint_variants.vcf.gz --missing-indv --out bcftools_variants/qc_missing
```
This creates `qc_missing.imiss` with columns:
- INDV — sample name
- N_DATA — sites with non-missing genotype calls
- N_GENOTYPES_FILTERED — sites filtered out for this sample
- N_MISS — missing genotype calls
- F_MISS — fraction missing (= N_MISS / (N_DATA + N_MISS))

Quick look at the table:
```bash
head bcftools_variants/qc_missing.imiss
```

Convert **fraction** to **percent** and list the worst samples (highest missing first):
```bash
# Tabular, percent missing, sorted descending
awk 'NR>1 {printf "%s\t%.2f\n", $1, $5*100}' bcftools_variants/qc_missing.imiss | sort -k2,2nr | column -t | head
```

Optional: compute the **mean missingness** across samples:
```bash
awk 'NR>1 {sum+=$5; n++} END {printf "Mean missing (%%): %.2f\n", (sum/n)*100}' bcftools_variants/qc_missing.imiss
```

**❓Questions (Missingness)**  
- Which sample has the highest missing rate? What is its percent missing?
- What threshold would you choose to exclude samples? (e.g., >10–20% missing)
- Are there samples that are outliers compared to the cohort mean?
- If one or two samples are very poor, what could be the causes (low coverage, library issues, contamination, wrong reference)?


### Section 7: Variant Filtering (with VCFtools)

After variant calling, the VCF file may contain many raw variants, including sequencing errors, rare artifacts, or positions with too much missing data.  
Filtering is a key step to obtain a **clean, high-confidence set of variants**.


#### Common Filters Explained

- **MAF (Minor Allele Frequency)**  
  - Definition: the frequency of the *less common allele* at a site.  
  - Example: If in 100 samples, 95 have allele A and 5 have allele G → MAF = 0.05 (5%).  
  - Why filter? Variants with very low MAF are often sequencing errors or not informative for population analyses.  
  - Typical threshold: **MAF ≥ 0.01 (1%)** or **MAF ≥ 0.05 (5%)**, depending on the project.  

- **Biallelic sites only**  
  - Many analyses assume only **two alleles** (reference and alternate).  
  - Multiallelic sites (more than 2 alleles) can complicate downstream analyses.  
  - Filter: `--min-alleles 2 --max-alleles 2`.

- **Missing data filter**  
  - If too many samples have missing genotypes at a variant, it may not be reliable.  
  - Example: `--max-missing 0.8` keeps only sites genotyped in at least 80% of samples.  

- **Remove indels**  
  - Indels (insertions/deletions) are harder to call reliably compared to SNPs.  
  - For training, we will focus on SNPs only.  


**7.1 Create a filtering job script**

Open a new file `vcf_filter.sh` with **vi**:

```bash
vi vcf_filter.sh
```

Press `i` to enter insert mode and paste:

```bash
#!/bin/bash
#SBATCH --job-name=vcf_filter
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH -o vcf_filter.out
#SBATCH -e vcf_filter.err

module load vcftools
module load bcftools

# Directories
VCF_DIR="/lisc/scratch/course/pgbiow/05_mapping_varriant_calling/bcftools_variants"
cd "$VCF_DIR" || exit 1

# Input joint VCF file
IN_VCF="joint_variants.vcf.gz"

# Apply filtering with VCFtools
vcftools --gzvcf "$IN_VCF" \
  --remove-indels \
  --min-alleles 2 --max-alleles 2 \
  --maf 0.01 \
  --max-missing 0.9 \
  --recode --recode-INFO-all \
  --out joint_variants.filtered

# Compress and index the filtered VCF
bgzip -f joint_variants.filtered.recode.vcf
bcftools index -t joint_variants.filtered.recode.vcf.gz

# Get new summary statistics
bcftools stats joint_variants.filtered.recode.vcf.gz > joint_variants.filtered.stats.txt

echo "Filtering complete. Results saved to:"
echo "$VCF_DIR/joint_variants.filtered.recode.vcf.gz"
echo "$VCF_DIR/joint_variants.filtered.stats.txt"
```

Save and exit (`ESC`, then `:wq`).

Submit the job:

```bash
sbatch vcf_filter.sh
```

**📍 Checkpoint:**  
After the job completes, you should have:
`joint_variants.filtered.recode.vcf.gz` → the filtered VCF file
`joint_variants.filtered.recode.vcf.gz.tbi` → index file
`filtered_stats.txt` → summary statistics after filtering

Please check the summary statistics after filtering:
```bash
head -n 30 bcftools_variants/filtered_stats.txt
```

**❓Questions**  
- How many SNPs were present before filtering, and how many remain after?
- Which filter (MAF, missing data, indels, biallelic) do you think removed the most variants?
- What problems might arise if we allow variants with too much missing data?
- Why is filtering out variants with very low MAF important for downstream analyses like PCA or GWAS?


## Conclusion

Congratulations! You have completed the genome mapping and variant calling workflow. You have learned how to take raw sequencing reads, align them to a reference genome, and identify genetic variants. You have also learned how to perform these tasks efficiently on an HPC cluster using Slurm and `for` loops.

This is a foundational workflow in genomics, and the skills you have learned today will be applicable to many different types of sequencing data analysis.

## You have completed the 1st Week! 

___
## Additional Tools for Mapping and Variant Calling

While in this training we used **BWA + Samtools + BCFtools**, other tools exist with different strengths, optimizations, and features.  
Here are some commonly used alternatives:

### Mapping Tools
- **Bowtie2** – very fast, suited for short reads and large genomes.  
  [Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)  

- **HISAT2** – fast and memory-efficient, often used for both DNA and RNA-seq alignment.  
  [HISAT2 Documentation](https://daehwankimlab.github.io/hisat2/)  

- **STAR** – very fast spliced aligner, primarily for RNA-seq data.  
  [STAR GitHub](https://github.com/alexdobin/STAR)  

- **Minimap2** – the go-to aligner for long reads (PacBio, Nanopore).  
  [Minimap2 Paper + GitHub](https://github.com/lh3/minimap2)  


### Variant Calling Tools
- **GATK (Genome Analysis Toolkit)** – industry standard for germline and somatic variant calling, with extensive best-practice pipelines.  
  [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us)  

- **FreeBayes** – haplotype-based variant detector, works well with pooled samples.  
  [FreeBayes GitHub](https://github.com/freebayes/freebayes)  

- **VarScan2** – robust for low-frequency variants and cancer studies.  
  [VarScan2 Manual](http://varscan.sourceforge.net/)  

- **DeepVariant** – Google’s deep learning-based variant caller, very accurate but computationally heavy.  
  [DeepVariant GitHub](https://github.com/google/deepvariant)  


### Useful Tutorials and Resources

- [Samtools Tutorial (Galaxy Training)](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/mapping/tutorial.html)  
- [A practical introduction to GATK 4](https://hpc.nih.gov/training/gatk_tutorial/)  
- [Variant calling with GATK](https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/)
- [NGS Analysis](https://learn.gencore.bio.nyu.edu/variant-calling/)
- [Genome mapping](https://genomics.sschmeier.com/ngs-mapping/))    
