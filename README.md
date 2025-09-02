# Verdant Bioinformatics Training 

This repository contains nine day-by-day pratical lessons for the Verdant bioinformatics training. Materials provide copy‑pasteable commands and Slurm job templates covering: Linux fundamentals; data QC; read mapping and variant calling; population structure and diversity; GWAS; ROH/linkage mapping; and RNA‑seq. All work is executed on LiSC HPC Cluster via the login node and Slurm. Access requires enabling the LiSC firewall (12‑hour window) and connecting via SSH to login02.lisc.univie.ac.at.

---

## Repository structure (Days 1–9)
The exact lesson files in this repository:

- **Day 1 — _Day 1. Introduction to Linux.md_**
- **Day 2 — _Day 2. Bioinformatics scripting essentials.md_**
- **Day 3 — _Day 3. Biological databases.md_**
- **Day 4 — _Day 4. Data QC and Experimental design.md_**
- **Day 5 — _Day 5. Mapping and variant calling.md_**
- **Day 6 — _Day 6. Diversity and Structure Analysis with SNPs.md_**
- **Day 7 — _Day 7. SNP applications in Breeding - GWAS.md_**
- **Day 8 — _Day 8. Runs of Homozygosity and Linkage Mapping.md_**
- **Day 9 — _Day 9. RNA-seq analysis.md_**

Each Day file contains copy‑pasteable commands and Slurm-ready snippets used in class.

---

## Daily access checklist (HPC - Cluster)
You’ll normally connect from outside Uni/MedUni networks. **Before every session/day** (or when your IP changes):

1. **Enable LiSC firewall access (12 h)**  
   Open: https://lisc.univie.ac.at/firewall/ → enter your **LiSC username** and confirm either via the **one-time link sent by e‑mail** or a **TOTP code** (recommended). Access lasts **12 hours**; renew as needed. Users on **dynamic IPs** should **terminate access after use** on the firewall page.

2. **SSH to the login node**  
   ```bash
   ssh <your_lisc_username>@login02.lisc.univie.ac.at
   ```

3. **Run heavy work via Slurm** (compute nodes only). Keep the login node for light tasks, editing, and submissions only.

---

## Slurm basics we use in the lessons
First, **discover and load modules** for the tools you need:

```bash
# discover what's available
module avail 2>&1 | head -n 30       # quick tree view
module spider bwa                     # search details & versions for a tool

# load specific versions (examples)
module load bwa/0.7.17 samtools/1.20 bcftools/1.20
module list                           # verify what's loaded

# sanity check the binaries
which bwa && bwa 2>&1 | head -n1
```

Then, **submit** batch jobs and **monitor progress**:

```bash
# submit
sbatch myjob.slurm

# watch queue
squeue -u $USER

# view accounting / finished jobs
sacct -j <jobid> --format=JobID,JobName%25,State,Elapsed,MaxRSS,AllocCPUS
```

A minimal template (`templates/`) will be referenced in the Day files; adapt `--partition`, `--account`, CPUs, memory, and time to your environment.

---

## Software used (by Day)

> Versions vary by cluster; anything in the listed major versions generally works. Items marked **(optional)** are alternates or advanced extras.

### Core HPC & CLI (all days)
- **Shell & Core**: `bash`, GNU coreutils, `grep`, `sed`, `awk`, `find`, `sort`, `uniq`, `xargs`, `tar`, `gzip`/`pigz`, `rsync`
- **Sessioning**: `screen` or `tmux`
- **Transfer/Networking**: `ssh`, `scp`, `sftp`, `wget`, `curl`
- **Version control**: `git`
- **Schedulers**: Slurm (`sbatch`, `srun`, `squeue`, `sacct`)
- **Modules**: Lmod/Environment Modules (`module`)

### Day 1–2 (Linux & scripting)
- Everything in **Core HPC & CLI**
- **Editors**: `nano`, `vim` (choose one)
- **Small utilities**: `jq` **(optional)**, `csvkit` **(optional)**, `parallel` **(optional)**
- **Python 3** (tiny helpers), **R** (quick plots)

### Day 3 (Biological databases)
- **NCBI SRA Toolkit**: `prefetch`, `fasterq-dump`
- **Entrez Direct (EDirect)**: `esearch`, `efetch`, `xtract`
- **MD5/SHA tools** for checksum validation
- **(optional)**: `datasets` CLI; **BLAST+** (quick sequence checks)

### Day 4 (QC & trimming)
- **fastqc**, **multiqc**, **fastp** (or **cutadapt** as alternative)
- **seqtk** **(optional)**

### Day 5 (Mapping & variant calling)
- **bwa** (short‑read alignment)
- **samtools** (BAM ops), **htslib** (bgzip/tabix)
- **picard** (duplicate marking)
- **bcftools** (pileup/call, merge, index)
- **vcftools** (summaries/filters)
- **(optional)**: `minimap2`, `bowtie2`

### Day 6 (Diversity & structure)
- **plink** (1.9/2.0) – LD pruning, PCA, QC
- **admixture** – ancestry components & CV
- **vcftools**, **bcftools** – filtering/conversions
- **R packages**: `tidyverse`, `data.table`, `ggplot2`, `cowplot`
- **(optional)**: **ANGSD**, **treemix**, **pixy**

### Day 7 (GWAS)
- **plink2** – association tests (linear/logistic), covariates
- **R packages**: `qqman`, `tidyverse`
- **(optional)**: **GEMMA**/**EMMAX** (mixed models), `plink --glm` (1.9)

### Day 8 (ROH & linkage mapping)
- **plink** – ROH (`--homozyg`)
- **(optional)**: **bcftools roh** (HMM‑based), **Lep‑MAP3**/**MSTmap** (advanced linkage mapping)

### Day 9 (RNA‑seq)
- **STAR** (spliced alignment) **or** **salmon** (quasi‑mapping)
- **subread**/**featureCounts** (read summarization)
- **R packages**: `DESeq2`, `edgeR`, `tximport`, `tidyverse`
- **(optional)**: `hisat2` + `stringtie`

### Helpful extras (any day)
- **bedtools** (interval ops), **bedops** **(optional)**
- **IGV** (desktop) for manual inspection
- **Apptainer/Singularity** **(optional)** for containerized runs

---

## Copying files between LiSC and your computer

Ensure your LiSC firewall access is active before transferring. Run these commands from your local machine (not from the cluster).

# Download a file (LiSC → local)
scp <your_lisc_username>@login02.lisc.univie.ac.at:/path/on/cluster/file.txt .

# Download a directory
scp -r <your_lisc_username>@login02.lisc.univie.ac.at:/path/on/cluster/RESULTS ./RESULTS

# Upload a file (local → LiSC)
scp ./localfile.txt <your_lisc_username>@login02.lisc.univie.ac.at:/path/on/cluster/

# Large/robust transfer with resume + progress
rsync -avhP <your_lisc_username>@login02.lisc.univie.ac.at:/path/on/cluster/RESULTS/ ./RESULTS/

# Interactive session
sftp <your_lisc_username>@login02.lisc.univie.ac.at

---
### Maintainer
Dr. Luiz Augusto Cauz dos Santos – University of Vienna 
