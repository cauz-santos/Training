# eDNA Microbiome Profiling with Kraken2 + Bracken + Pavian

This guide provides a **complete pipeline** for analyzing eDNA samples sequenced with GBS-like reads using **Kraken2** and **Bracken**, followed by visualization with **Pavian** in R.

---

## 1. Install Dependencies

We will use conda (Miniconda/Anaconda).

```bash
# Create environment
conda create -n kraken2_env -c bioconda -c conda-forge kraken2 bracken krona fastqc fastp -y

# Activate environment
conda activate kraken2_env
```

This installs:
- **Kraken2** (classifier)
- **Bracken** (abundance estimation)
- **Krona** (HTML visualization)
- **FastQC/fastp** (quality control)

---

## 2. Download the Database

The best option is the **PlusPF database**, which includes:
- Bacteria, Archaea, Viruses, Plasmids, Protozoa, Fungi, and Human (for host filtering).  
- Size: ~450 GB disk, needs ~200 GB RAM.

```bash
# Create database folder
mkdir -p kraken2_db/pluspf
cd kraken2_db/pluspf

# Download libraries
kraken2-build --download-library bacteria --db pluspf_db
kraken2-build --download-library archaea --db pluspf_db
kraken2-build --download-library viral --db pluspf_db
kraken2-build --download-library plasmid --db pluspf_db
kraken2-build --download-library protozoa --db pluspf_db
kraken2-build --download-library fungi --db pluspf_db
kraken2-build --download-library human --db pluspf_db

# Build database
kraken2-build --build --db pluspf_db

# Return to main directory
cd ../../
```

If resources are limited, consider using the **Standard DB (~100 GB)** or **MiniKraken (~8 GB)**.

---

## 3. Quality Control

Run QC and trimming on raw FASTQ files:

```bash
fastqc sample_R1.fastq sample_R2.fastq

fastp -i sample_R1.fastq -I sample_R2.fastq \
      -o clean_R1.fastq -O clean_R2.fastq \
      -h fastp.html -j fastp.json
```

---

## 4. Classify Reads with Kraken2

```bash
kraken2 --db kraken2_db/pluspf/pluspf_db \
  --threads 16 \
  --paired clean_R1.fastq clean_R2.fastq \
  --report sample.kraken.report \
  --output sample.kraken.out
```

- `sample.kraken.report`: summary table of taxa  
- `sample.kraken.out`: classification for each read  

---

## 5. Refine Abundance with Bracken

Bracken refines counts into better abundance estimates:

```bash
bracken -d kraken2_db/pluspf/pluspf_db \
  -i sample.kraken.report \
  -o sample.bracken.species \
  -r 150 -l S
```

- `-r 150`: read length (adjust if yours is different)  
- `-l S`: taxonomic level (S = species, G = genus, etc.)  

---

## 6. Visualization with Krona

```bash
ktImportTaxonomy sample.bracken.species -o sample.krona.html
```

Open `sample.krona.html` in your browser.

---

## 7. Visualization with Pavian (R)

Install and run Pavian:

```R
# In R
install.packages("remotes")
remotes::install_github("fbreitwieser/pavian")

library(pavian)
pavian::runApp(port=5000)
```

Then, open your browser at [http://127.0.0.1:5000](http://127.0.0.1:5000) and upload:  
- `sample.kraken.report`  
- or `sample.bracken.species`  

You can explore interactive barplots, heatmaps, and tables.
