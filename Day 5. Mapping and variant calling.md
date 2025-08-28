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

We will be working with a set of up to 10 single-end FASTQ files and a reference genome. The workflow is as follows:

1.  **Setup:** Prepare the environment and data.
2.  **Reference Genome Indexing:** Create an index for the reference genome to enable fast alignment.
3.  **Read Mapping:** Align the FASTQ reads to the indexed reference genome.
4.  **Alignment Post-processing:** Convert SAM to BAM, sort, and index the alignment files.
5.  **Quality Control of Alignments:** Check mapping rates and coverage.
6.  **Variant Calling:** Identify genetic variants from the aligned reads.
7.  **Variant Filtering and Annotation:** Filter raw variants to improve accuracy and interpret the results.

### The Data

*   **Reference Genome of *Elaeis guineensis*:** `GCF_000442705.2_EG11_genomic.fna` (Size: 1.8 GB)
*   **Input Reads:** 12 single-end FASTQ files (e.g., `DRR070477.fastq.gz`, `DRR070483.fastq.gz`, `DRR070491.fastq.gz`, `DRR070494.fastq.gz`, `DRR070497.fastq.gz`, `SRR14510793.fastq.gz`, `SRR14510807.fastq.gz`, `DRR070482.fastq.gz`, `DRR070488.fastq.gz`, `DRR070492.fastq.gz`, `DRR070496.fastq.gz`, `DRR070498.fastq.gz`, `SRR14510795.fastq.gz`)

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

    ```bash
    mkdir -p day5_mapping_variant_calling
    cd day5_mapping_variant_calling
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

Now, let's create the `index_genome.sh` file:

```bash
#!/bin/bash
#SBATCH --job-name=bwa_index
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -o bwa_index.out
#SBATCH -e bwa_index.err

# Load the BWA module
module load bwa

# Define the reference genome file
GENOME="GCF_000442705.2_EG11_genomic.fna"

# Run BWA index
echo "Indexing the reference genome: $GENOME"
bwa index $GENOME

echo "Indexing complete."
```

**Submit the job:**

```bash
sbatch index_genome.sh
```

**❓ Question:** What do the `#SBATCH` directives in the script mean? Why are they important for running jobs on a cluster?

Once the job is complete, you will see several new files in your directory with extensions like `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`. These are the BWA index files.



### Section 3: Read Mapping and Post-processing

Now that we have our indexed genome, we can align our FASTQ reads. We will create a single Slurm script that uses a `for` loop to process all our FASTQ files.

**Create a file named `map_reads.sh` using `vi`:**

```bash
#!/bin/bash
#SBATCH --job-name=read_mapping
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o read_mapping.out
#SBATCH -e read_mapping.err

# Load necessary modules
module load bwa
module load samtools

# Define the reference genome
GENOME="GCF_000442705.2_EG11_genomic.fna"

# Loop through all FASTQ files in the current directory
for fq in *.fastq.gz
do
    # Get the base name of the file for output naming
    base=$(basename "$fq" .fastq.gz)

    echo "Processing sample: $base"

    # 1. Align reads with BWA-MEM
    # The output is piped directly to samtools to convert to BAM
    echo "  Aligning reads..."
    bwa mem -t 8 $GENOME "$fq" | samtools view -bS - > "${base}.bam"

    # 2. Sort the BAM file
    # Sorting is required for many downstream tools, including variant calling
    echo "  Sorting BAM file..."
    samtools sort "${base}.bam" -o "${base}.sorted.bam"

    # 3. Index the sorted BAM file
    # Indexing allows for fast random access to the alignment data
    echo "  Indexing sorted BAM file..."
    samtools index "${base}.sorted.bam"

    # Clean up intermediate files
    rm "${base}.bam"

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

You can do this in a `for` loop on the command line:

```bash
for bam in *.sorted.bam
do
    echo "Statistics for $bam:"
    samtools flagstat "$bam"
    echo "--------------------"
done
```

**Look for the line that says `... + ... mapped (...%)`.** This is your overall mapping rate.

**❓ Question:** What would you consider a good mapping rate? What might be the cause of a low mapping rate?


**4.2: Checking Coverage**  

Coverage (or read depth) is the number of reads that align to, or "cover," a specific position in the genome. Adequate coverage is essential for confident variant calling. Low coverage can lead to false negatives (missing real variants), while extremely high coverage might indicate PCR duplicates or repetitive regions.

We can use `samtools depth` or `samtools coverage` to assess coverage.

**Calculate average coverage across the genome:**

```bash
# This command calculates the depth at each position and then uses awk to average it.
# This can be slow for a large genome.
samtools depth sample1.sorted.bam | awk 
'{sum+=$3} END { print "Average coverage = ",sum/NR}
'
```

**A more comprehensive view with `samtools coverage`:**

This command provides a summary of coverage statistics, including the percentage of the genome covered at different depths.

```bash
# Run on a single BAM file
samtools coverage sample1.sorted.bam
```

**❓ Question:** What is the difference between breadth of coverage and depth of coverage? Why are both important?

**4.3: Inspecting Alignments Visually**

Sometimes, it's helpful to look at the alignments in a specific region, especially if you have a gene of interest. `samtools tview` provides a text-based alignment viewer directly in your terminal.

```bash
# You need the reference genome file here
samtools tview sample1.sorted.bam GCF_000442705.2_EG11_genomic.fna
```

Navigate to a specific region by pressing `g` and entering the coordinates (e.g., `chromosome_name:10000`).

**❓ Question:** What do the different characters and colors in the `tview` output represent?



### Section 5: Variant Calling

Now that we have high-quality alignments, we can proceed to variant calling. This process identifies positions where the aligned reads consistently differ from the reference genome.

We will use `bcftools mpileup` and `bcftools call`.

*   `bcftools mpileup`: Summarizes the base calls of aligned reads at each position in the genome.
*   `bcftools call`: Calls SNPs and indels from the `mpileup` output.

**Create a Slurm script named `variant_calling.sh` using `vi`:**

```bash
#!/bin/bash
#SBATCH --job-name=variant_calling
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH -o variant_calling.out
#SBATCH -e variant_calling.err

module load bcftools

GENOME="GCF_000442705.2_EG11_genomic.fna"

# Create list of BAM files
ls *.sorted.bam > bamlist.txt

# Joint calling
bcftools mpileup -Ou -f $GENOME -b bamlist.txt | \
bcftools call -mv -Oz -o joint_variants.vcf.gz

# Index VCF
bcftools index joint_variants.vcf.gz

echo "Variant calling complete."
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
bcftools view joint_variants.vcf.gz | head -n 50
```

**The VCF file has two main parts:**

*   **Header:** Lines starting with `##` describe the file format, the reference genome, and the annotations used in the file (INFO, FORMAT fields).
*   **Variant entries:** Each line represents a variant, with columns for chromosome, position, reference allele, alternate allele, quality, and sample-specific information.

**❓ Question:** What is the difference between the `INFO` and `FORMAT` fields in a VCF file?

**6.2: VCF Statistics**

`bcftools stats` provides a wealth of information about your VCF file, including the number of SNPs and indels, transition/transversion ratio, and more.

```bash
bcftools stats joint_variants.vcf.gz > joint_variants.vcf.stats
```

You can then view the generated statistics file:

```bash
less joint_variants.vcf.stats
```

**❓ Question:** What is the transition/transversion (Ts/Tv) ratio, and why is it a useful quality metric for SNP calls?

**6.3: Checking Missing Data (per sample)**

Missing genotypes can bias downstream population analyses. Let’s quantify **percent missing per sample** on the *joint VCF*.

We’ll use **VCFtools**

```bash
#First load the module
module load vcftools

# Compute per-individual missingness
vcftools --gzvcf joint_variants.vcf.gz --missing-indv --out qc_missing
```
This creates `qc_missing.imiss` with columns:
- INDV — sample name
- N_DATA — sites with non-missing genotype calls
- N_GENOTYPES_FILTERED — sites filtered out for this sample
- N_MISS — missing genotype calls
- F_MISS — fraction missing (= N_MISS / (N_DATA + N_MISS))

Quick look at the table:
```bash
head qc_missing.imiss
```

Convert **fraction** to **percent** and list the worst samples (highest missing first):
```bash
# Tabular, percent missing, sorted descending
awk 'NR>1 {printf "%s\t%.2f\n", $1, $5*100}' qc_missing.imiss | sort -k2,2nr | column -t | head
```

Optional: compute the **mean missingness** across samples:
```bash
awk 'NR>1 {sum+=$5; n++} END {printf "Mean missing (%%): %.2f\n", (sum/n)*100}' qc_missing.imiss
```

**❓Questions (Missingness)**  
- Which sample has the highest missing rate? What is its percent missing?
- What threshold would you choose to exclude samples? (e.g., >10–20% missing)
- Are there samples that are outliers compared to the cohort mean?
- If one or two samples are very poor, what could be the causes (low coverage, library issues, contamination, wrong reference)?

### Section 7: Variant Filtering

Raw variant calls often contain false positives due to sequencing errors, mapping errors, or other artifacts. Filtering is a critical step to increase the confidence in your variant set.

We will use `bcftools filter` to apply some common filters.

**Example of filtering:**

Let's filter for variants with a quality score greater than 30 and a read depth of at least 10.

```bash
bcftools filter -i 'QUAL > 30 && DP > 10' -Oz -o sample1.filtered.vcf.gz sample1.vcf.gz
```

**Explanation of the filter expression:**

*   `-i '...'`: Includes sites for which the expression is true.
*   `QUAL > 30`: Selects variants with a quality score above 30.
*   `DP > 10`: Selects variants where the read depth is greater than 10.

**❓ Question:** What other filtering criteria might be useful? Think about strand bias, mapping quality, and other fields available in the VCF file.

**Advanced Filtering:**

You can create more complex filtering rules. For example, to filter based on attributes in the `INFO` field:

```bash
# Example: Filter based on mapping quality bias (MQ) and strand bias (SB)
bcftools filter -e 'MQ < 40 || SB > 0.1' -Oz -o sample1.advanced_filtered.vcf.gz sample1.filtered.vcf.gz
```

*   `-e '...'`: Excludes sites for which the expression is true.

---

## Conclusion

Congratulations! You have completed the genome mapping and variant calling workflow. You have learned how to take raw sequencing reads, align them to a reference genome, and identify genetic variants. You have also learned how to perform these tasks efficiently on an HPC cluster using Slurm and `for` loops.

This is a foundational workflow in genomics, and the skills you have learned today will be applicable to many different types of sequencing data analysis.
