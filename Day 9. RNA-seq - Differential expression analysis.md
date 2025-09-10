## Week 2: Applied Bioinformatics for Genomics and Breeding

## Day 9: RNA-seq Data Analysis

Day 9 is dedicated to RNA sequencing (RNA-seq) data analysis, a powerful high-throughput technology used to measure gene expression levels, discover novel transcripts, and identify gene fusions. We will cover the entire workflow, from experimental design and pre-processing to mapping, quantification, and differential expression analysis.

**RNA-seq Experimental Design and Pre-processing**

As discussed on Day 5, robust experimental design is critical for RNA-seq studies. This session will reiterate key design principles specific to RNA-seq and then dive into the initial steps of data pre-processing.

**Key Considerations for RNA-seq Experimental Design:**

*   **Biological Replicates:** Absolutely essential for statistical power and distinguishing true biological variation from technical noise. Aim for at least 3-5 biological replicates per condition.
*   **Sequencing Depth:** The number of reads per sample. Deeper sequencing allows for better detection of lowly expressed genes and more accurate quantification. Typical depths range from 10-50 million reads per sample, depending on the complexity of the transcriptome and research question.
*   **Read Length:** Longer reads can improve alignment accuracy, especially in regions with repetitive sequences or alternative splicing.
*   **Paired-end vs. Single-end:** Paired-end reads provide more information for alignment and transcript assembly, especially for identifying splice junctions.
*   **RNA Quality and Integrity:** High-quality RNA is crucial. RNA Integrity Number (RIN) scores are commonly used to assess RNA degradation. Degraded RNA can lead to biased results.
*   **Batch Effects:** As mentioned, samples processed at different times or by different personnel can introduce systematic biases. Randomization of samples across batches is important.


## Learning outcomes
By the end of this session you will be able to:
- Build a **STAR** genome index and align paired‑end RNA‑seq reads (demo sample).
- Generate a **gene‑level count matrix** using **featureCounts** (strandedness‑aware).
- Perform **differential expression (DE)** with **edgeR via Trinity** helper scripts.
- Run **GO enrichment** (GOseq via Trinity or g:Profiler) and interpret pathways.
- Tie RNA‑seq findings to **plant breeding** decisions (candidate pathways/genes for P‑use efficiency).

---

## Dataset used in this class (Oil palm, Pi‑starvation, 28 days)
**Species/tissue:** *Elaeis guineensis* (oil palm), roots.  
**Design:** 6 libraries at **28‑day** timepoint — **Pi‑starved (−P)** vs **Pi‑sufficient (+P)**; **3 biological replicates per group**.  
**Source:** Kong *et al.* (2021) “Comparative transcriptome analysis reveals novel insights into transcriptional responses to phosphorus starvation in oil palm root.” (BioProject PRJNA673667).  
**Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7863428/

| Group | Time | Tissue | Rep | SRA Run | Sample name |
|---|---|---|---:|---|---|
| −P (Pi‑starved) | 28 d | root | 1 | SRR12959886 | Pminus_28d_rep1 |
| −P (Pi‑starved) | 28 d | root | 2 | SRR12959894 | Pminus_28d_rep2 |
| −P (Pi‑starved) | 28 d | root | 3 | SRR12959895 | Pminus_28d_rep3 |
| +P (Pi‑sufficient) | 28 d | root | 1 | SRR12959887 | Pplus_28d_rep1 |
| +P (Pi‑sufficient) | 28 d | root | 2 | SRR12959888 | Pplus_28d_rep2 |
| +P (Pi‑sufficient) | 28 d | root | 3 | SRR12959889 | Pplus_28d_rep3 |

> **Note:** We **start from trimmed FASTQs** and **pre‑mapped BAMs** (except one demo sample we’ll map live).


**`metadata/rnaseq_samples.txt` (Trinity format)**
```
Pminus    Pminus_28d_rep1,Pminus_28d_rep2,Pminus_28d_rep3
Pplus     Pplus_28d_rep1,Pplus_28d_rep2,Pplus_28d_rep3
```

**`metadata/rnaseq_sample_sheet.csv`**
```csv
sample_id,group,batch
Pminus_28d_rep1,minus,B1
Pminus_28d_rep2,minus,B1
Pminus_28d_rep3,minus,B1
Pplus_28d_rep1,plus,B1
Pplus_28d_rep2,plus,B1
Pplus_28d_rep3,plus,B2
```

> **Breeding relevance:** Phosphorus is costly and often limiting; identifying transcriptional responses underpinning **P‑use efficiency (PUE)** enables **marker development**, **line selection**, and **reduced fertilizer** strategies.

---

## Cluster modules (LISC)
```bash
module purge
module load Trinity/2.15.2-foss-2023a   # or: module load trinityrnaseq/2.15.2-5.40.2-5.3.28
module load STAR/2.7.10a                # use the closest STAR available
module load SAMtools/1.20
module load Subread/2.0.3               # provides featureCounts
export TRINITY_HOME=${TRINITY_HOME:-$EBROOTTRINITY}
```

---


First create a folder for the file outputs of day 9 in your home directory:
   ```bash
   mkdir 09_rnaseq_expression
   ```
Now, enter the folder with the command `cd`

Then create some subfolders, for each specific analysis we will perform:
   ```bash
   mkdir reference star trinity go_terms
   ```

### 1) STAR genome index (Slurm)

**What are we doing and why?**  
We will **build a STAR genome index** from a reference **FASTA** and **GTF**. STAR needs this index to align reads quickly and to detect **splice junctions** accurately.  
>**Breeding relevance:** robust splice-aware alignment preserves signals from genes and isoforms involved in **phosphorus-use efficiency (PUE)** and other stress responses—exactly the biology we want to prioritize for selection, marker design, and follow-up validation.

**Create the script with `vi`:**  
Open a new file:
```bash
vi 01_star_index.sh
```

When `vi` opens, press `i` to enter insert mode.

Paste the script below :

```bash
#!/usr/bin/env bash
#SBATCH -p standard
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH -t 01:00:00
#SBATCH -J day9_star_index
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail

# --- Load modules ---
module purge
module load Trinity/2.15.2-foss-2023a
module load STAR

# --- Define paths ---
GENOME_DIR="/lisc/scratch/course/pgbiow/data/genomes"
WORK_REF="reference"

# --- Copy reference genome & annotations into your working folder ---
mkdir -p "$WORK_REF"
cp "${GENOME_DIR}/Elaeis_guineensis_genomic.fna" "$WORK_REF/reference.fa"
cp "${GENOME_DIR}/Elaeis_guineensis_genomic.gtf" "$WORK_REF/annotation.gtf"
# (optional: keep gff too if you need)
cp "${GENOME_DIR}/Elaeis_guineensis_genomic.gff" "$WORK_REF/annotation.gff"

# --- Define input files for STAR ---
REF="$WORK_REF/reference.fa"
GTF="$WORK_REF/annotation.gtf"
IDX="$WORK_REF/STAR_index"

READLEN=150          # set to your read length (100 or 150)
SJDB=$((READLEN-1))

# --- Make index directory ---
mkdir -p "$IDX"

# --- Run STAR genomeGenerate ---
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir "$IDX" \
     --genomeFastaFiles "$REF" \
     --sjdbGTFfile "$GTF" \
     --sjdbOverhang $SJDB

echo "STAR index built successfully in $IDX"
```
To save and exit `vi`: press `Esc`, then type `:wq` and press `Enter`.

Submit:
```bash
sbatch 01_star_index.sh
```

**How to check success:**  
STAR will populate ref/STAR_index/ with many small files (e.g., SA, Genome, sjdbList.out.tab).

Inspect the Slurm log for any error messages:
```bash
tail -n +1 logs/day9_star_index_*.out
tail -n +1 logs/day9_star_index_*.err
```

> *Breeding note:* Accurate splice‑aware mapping preserves signals from **splicing/TF regulation** linked to PUE traits.


### 2) Map **one** demo sample (Slurm)
We will align **one** paired-end RNA-seq sample to the STAR index you just built. This produces a **coordinate-sorted, indexed BAM** plus STAR logs. Mapping only one sample live keeps the session fast while still showing you every critical decision (flags, resources, and checks).  
>**Breeding relevance:** consistent, high-quality alignment across treatments (−P vs +P) is essential to avoid false differential expression and to capture splice-aware signals in genes underlying **phosphorus-use efficiency (PUE)**.

#### What to look for in the output
- `Aligned.sortedByCoord.out.bam` (+ `.bai` index) for the demo sample  
- `Log.final.out` with key metrics:  
  - **Number of input reads**  
  - **Uniquely mapped %** (aim high)  
  - **% mapped to multiple loci** (too high can indicate repeats/contamination)  
  - **Mismatch rate per base** (quality/adapter issues)  
  - **Chimeric reads** (usually low for bulk RNA-seq)

**Create the script with `vi`:**  
1) Open the file:
```bash
vi 02_star_map_one.sh
```
2) **Press `i`** to enter *insert* mode.  
3) Paste the script.  
4) Save & quit: **`Esc`**, then **`:wq`**.

```bash
#!/usr/bin/env bash
#SBATCH -p standard
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 00:35:00
#SBATCH -J day9_star_map1
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail

module purge
module load Trinity/2.15.2-foss-2023a
module load STAR/2.7.10a
module load SAMtools/1.20

# Paths
IDX="/lisc/scratch/course/pgbiow/09_rnaseq_expression/reference/STAR_index"
FASTQ_DIR="/lisc/scratch/course/pgbiow/data/RNAseq"
OUT_BASE="/lisc/scratch/course/pgbiow/09_rnaseq_expression/star"

# Sample (change if needed)
SAMPLE=${SAMPLE:-Pminus_28d_rep1}
R1=${FASTQ_DIR}/${SAMPLE}_R1.paired.fastq
R2=${FASTQ_DIR}/${SAMPLE}_R2.paired.fastq
OUT=${OUT_BASE}/${SAMPLE}
mkdir -p "$OUT"

# Run STAR
STAR --runThreadN 16 \
  --genomeDir "$IDX" \
  --readFilesIn "$R1" "$R2" \
  --outFileNamePrefix "$OUT/" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS nM XS \
  --twopassMode Basic

# Index BAM
samtools index -@ 8 ${OUT}/Aligned.sortedByCoord.out.bam

# Quick mapping summary
echo "==== ${SAMPLE} mapping summary ===="
grep -E "Number of input reads|Uniquely mapped|% of reads mapped to multiple|Mismatch rate" ${OUT}/Log.final.out || true
```

Submit:
```bash
SAMPLE=Pminus_28d_rep1 sbatch 02_star_map_one.sh
```

After submission, we’ll **inspect `Log.final.out`** together and discuss whether the metrics support reliable downstream DE analysis.
*(All other samples are pre‑mapped and already under `star/`.)*


### 3) featureCounts → counts matrix (Slurm)

### What are we doing and why?
We will convert the per-sample **aligned BAM files** (from STAR) into a single **gene-level count matrix** using **featureCounts** (part of Subread). This step assigns read pairs to **exons** and sums them by **gene_id** from the GTF.  
>**Breeding relevance:** accurate gene counts are the backbone of **edgeR** differential expression. This is where nutrient/stress pathways (e.g., phosphate transport, root development) first become **quantifiable** for ranking candidates and informing selection decisions.

### Inputs and outputs
- **Inputs**
  - Sorted BAMs for all samples: `star/*/Aligned.sortedByCoord.out.bam`
  - Gene annotation (GTF): `ref/annotation.gtf`
- **Outputs**
  - Raw featureCounts table: `counts/featureCounts.txt`
  - Clean **genes × samples** matrix (TSV): `counts/counts_matrix.tsv`
  - `gene_lengths.txt` (for GOseq length bias correction): `counts/gene_lengths.txt`

### About strandedness (`-s`)
- `-s 0` = **unstranded**, `-s 1` = **forward**, `-s 2` = **reverse** (common for dUTP/Illumina TruSeq Stranded).  
- Using the wrong value can **halve** your counts or inflate antisense signal. If unsure, confirm with your library kit notes or a quick `RSeQC infer_experiment.py` test.  
- In this course, set it via an environment variable at submission (see **Submit** section).

**Create the script with `vi`:**  
1) Open the file:
```bash
vi 03_featurecounts.sh
```
2) **Press `i`** to enter *insert* mode.  
3) Paste the script below.  
4) Save & quit: **`Esc`**, **`:wq`**.

```bash
#!/usr/bin/env bash
#SBATCH -p standard
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 00:25:00
#SBATCH -J day9_featureCounts
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail

module purge
module load Trinity/2.15.2-foss-2023a
module load subread

# Reference annotation
GTF="/lisc/scratch/course/pgbiow/09_rnaseq_expression/reference/annotation.gtf"

# Output directory
OUT="/lisc/scratch/course/pgbiow/09_rnaseq_expression/counts"
mkdir -p "$OUT"

# Collect all BAMs
ls /lisc/scratch/course/pgbiow/data/RNAseq/Mapped/*.bam > ${OUT}/bam_list.txt

# Set strandedness: 0=unstranded, 1=forward, 2=reverse
STRAND=${STRAND:-0}

featureCounts -T 16 -p -B -C -s ${STRAND} \
  -t exon -g gene_id \
  -a "$GTF" \
  -o ${OUT}/featureCounts.txt \
  $(cat ${OUT}/bam_list.txt)

# Clean genes x samples matrix
python - << 'PY'
import pandas as pd
fc = pd.read_csv('/lisc/scratch/course/pgbiow/09_rnaseq_expression/counts/featureCounts.txt', sep='\t', comment='#')
cts = fc.iloc[:, [0] + list(range(6, fc.shape[1]))].copy()
cts.rename(columns={'Geneid':'gene_id'}, inplace=True)
cts.to_csv('/lisc/scratch/course/pgbiow/09_rnaseq_expression/counts/counts_matrix.tsv', sep='\t', index=False)
# gene lengths for GOseq
fc[['Geneid','Length']].rename(columns={'Geneid':'gene_id','Length':'length'}) \
  .to_csv('/lisc/scratch/course/pgbiow/09_rnaseq_expression/counts/gene_lengths.txt', sep='\t', index=False, header=False)
PY
```

Submit:
```bash
# Example: reverse-stranded libraries
STRAND=2 sbatch 03_featurecounts.sh
```

To clean the sample names, removing unnecessary information, we can use:
```bash
awk 'BEGIN{OFS=FS="\t"} NR==1{for(i=2;i<=NF;i++) {sub(".*/", "", $i); sub(/_Aligned\.sortedByCoord\.out\.bam$/, "", $i)}; print} NR>1{print}' counts_matrix.tsv > counts_matrix_clean.tsv
```

> *Breeding note:* Gene‑level counts quantify **stress and nutrient pathways** driving resilience and yield stability under low‑P.

### After featureCounts: how to inspect your outputs

#### 1) Inspect the summary table
The file `featureCounts.txt.summary` reports, **per sample**, how many fragments were `Assigned` (counted to a gene) and how many were `Unassigned_*` for various reasons.

```bash
# Quick view
less -S counts/featureCounts.txt.summary
```
**What to look for:**
- `Assigned` should be the **largest** row across samples.
- Large `Unassigned_NoFeatures` → GTF mismatch (wrong genome/annotation).
- Large `Unassigned_Strand` → wrong `-s` value (strandedness).
- Large `Unassigned_MultiMapping` → repetitive regions or liberal aligner settings.

**Compute Assigned % per sample (quick check):**
```bash
awk 'NR==1{for(i=2;i<=NF;i++) S[i]=$i; next} {for(i=2;i<=NF;i++) T[i]+=$i} $1=="Assigned"{for(i=2;i<=NF;i++) A[i]=$i}
     END{for(i=2;i<=NF;i++){pct=(A[i]/T[i])*100; printf "%s\t%.2f%% Assigned\n", S[i], pct}}' counts/featureCounts.txt.summary
```

#### 2) Spot-check the clean counts matrix
The file `counts/counts_matrix.tsv` is **genes × samples** and is what you’ll feed into Trinity/edgeR.

```bash
# Header line (gene_id and sample names, check order!)
head -n 1 counts/counts_matrix_clean.tsv

# First 10 data rows
head counts/counts_matrix_clean.tsv | column -t
```

**Check total library sizes (sum of counts) per sample:**
```bash
awk -F'\t' 'NR==1{for(i=2;i<=NF;i++) name[i]=$i; next}
            {for(i=2;i<=NF;i++) lib[i]+=$i}
            END{for(i=2;i<=NF;i++) printf "%s\t%d\n", name[i], lib[i]}' counts/counts_matrix_clean.tsv | column -t
```

**Sanity check for all-zero genes (should be few or none):**
```bash
awk -F'\t' 'NR>1{sum=0; for(i=2;i<=NF;i++) sum+=$i; if(sum==0) print $1}' counts/counts_matrix_clean.tsv | head
```


**Common red flags & fixes**
- **Assigned < 60%** across samples → check `-s` and GTF/genome match.
- **One sample much smaller library size** → expect more dispersion; check `Log.final.out` and earlier QC.

### 4) Differential expression with **edgeR via Trinity** (Slurm)
#### What are we doing and why?
We’ll use Trinity’s wrappers for **edgeR** to run a clean two-group DE analysis (`Pminus` vs `Pplus`) with **TMM normalization**, generate per-contrast DE tables, and optional QC plots.  
**Breeding relevance:** DE highlights **candidate genes/TFs** and **pathways** that could be leveraged for **marker development**, **introgression**, or **editing** to improve PUE.

**Prepare `samples.txt` (Trinity format)**  
Create the file with `vi`:
```bash
vi trinity/samples.txt
```
Press **`i`**, paste, then **`Esc` → `:wq`**:
```
Pminus	Pminus_28d_rep1
Pminus	Pminus_28d_rep2
Pminus	Pminus_28d_rep3
Pplus	Pplus_28d_rep1
Pplus	Pplus_28d_rep2
Pplus	Pplus_28d_rep3
```
**Create the DE script with `vi`**  
```bash
vi 04_trinity_edger.sh
```
Press **`i`**, paste, then **`Esc` → `:wq`**:

```bash
#!/usr/bin/env bash
#SBATCH -p standard
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 00:45:00
#SBATCH -J day9_edgeR
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

# Load required software
module purge
module load Trinity/2.15.2-foss-2023a
module load R  # R must include edgeR, gplots, pheatmap, etc.

# Set Trinity differential expression directory
export TRINITY_HOME=${TRINITY_HOME:-$EBROOTTRINITY}
TRINITY_DE="${TRINITY_HOME}/Analysis/DifferentialExpression"

# Define output folder
OUTDIR="trinity"
mkdir -p "$OUTDIR"

# Print session info for reproducibility
echo "===== R session info ====="
Rscript -e "sessionInfo(); library(edgeR)" || { echo "ERROR: edgeR not available."; exit 1; }

echo "===== Running edgeR DE analysis ====="
$TRINITY_DE/run_DE_analysis.pl \
  --matrix counts/counts_matrix_clean.tsv \
  --method edgeR \
  --samples_file trinity/samples.txt \
  --output "$OUTDIR" \
  --min_reps_min_cpm 2,1

# Run clustering, PCA, and heatmaps on the DE results
echo "===== Running clustering and PCA (FDR ≤ 0.05, |log2FC| ≥ 1) ====="
cd "$OUTDIR"

$TRINITY_DE/analyze_diff_expr.pl \
  --matrix ../counts/counts_matrix_clean.tsv \
  --samples ../trinity/samples.txt \
  -P 0.05 \
  -C 1 \
  --output diffexpr_plots

echo "===== All done! Check 'trinity/diffexpr_plots/' for heatmaps and PCA ====="
```
Submit:
```bash
sbatch 04_trinity_edger.sh
```

**Outputs of interest**
- `edger_trinity/edgeR.DE_results/*DE_results` (per-contrast gene tables with logFC, FDR)
- `edger_trinity/diffexpr_plots/` (PCA, heatmaps, volcano/MA)

> *Breeding note:* DE highlights candidate genes/TFs for **marker design**, **introgression**, or **editing** to improve PUE.

### Outputs to Inspect — Trinity edgeR Differential Expression

This is how to quickly sanity-check what Trinity/edgeR produced and grab the numbers you need.

**1) Find your contrast name(s)**  
Trinity names outputs by contrast (e.g., `Pminus_vs_Pplus`). List them:
```bash
ls trinity/*DE_results | sed 's#.*/##' | sed 's/\.DE_results.*$//' | sort -u
```
> The **first** group in `A_vs_B` is the one considered **“up”** in `A_vs_B.DE_up.genes`.


**2) Inspect DE result tables**    
Each contrast has a `...DE_results` table with columns like `logFC`, `logCPM`, `PValue`, `FDR`.

***Preview a table (replace CONTR):***    
```bash
CONTR="counts_matrix_clean.tsv.Pminus_vs_Pplus"
DETAB="trinity/${CONTR}.edgeR.DE_results"
head -n 20 "$DETAB" | column -t
```

***Top 20 by lowest FDR (header-aware):***  
```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="FDR") f=i; print; next}
            NR>1{print $f "\t" $0}' "$DETAB" \
| sort -t$'\t' -k1,1g \
| cut -f2- \
| head -n 20 \
| column -t
```

***Top 20 by largest |logFC| (header-aware):***  
```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="logFC") l=i; print; next}
            NR>1{lf=$l; if(lf<0) lf=-lf; print lf "\t" $0}' "$DETAB" \
| sort -t$'\t' -k1,1gr \
| cut -f2- \
| head -n 20 \
| column -t
```

***Count significant DE genes at FDR ≤ 0.05:***  
```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="FDR") f=i; next} NR>1 && $f!~"NA" && $f<=0.05{c++} END{print c+0}' "$DETAB"
```


**3) Quick look at “up”/“down” gene lists**  
Trinity also writes simple lists of gene IDs:
```bash
# All DEGs with FDR ≤ 0.05 and |log2FC| ≥ 1
less trinity/${CONTR}.edgeR.DE_results.P0.05_C1.DE.subset

# Upregulated in Pminus
less trinity/${CONTR}.edgeR.DE_results.P0.05_C1.Pminus-UP.subset

# Upregulated in Pplus
less trinity/${CONTR}.edgeR.DE_results.P0.05_C1.Pplus-UP.subset
```
> **Interpretation:** `A_vs_B.DE_up.genes` are higher in **A** than **B**; `DE_down.genes` are lower in **A** than **B**.


**5) View plots (PCA, heatmaps, MA, Volcano)**  
List the files:
```bash
ls edger_trinity/diffexpr_plots/
```
If you’re on a headless cluster, copy them to your laptop:
```bash
# From your laptop/desktop terminal
scp -r USER@login.cluster:/work/USER/day09_rnaseq/edger_trinity/diffexpr_plots .
```
Then open the PDFs/PNGs locally.

**What to look for:**  
- **PCA**: replicates from the same group should cluster together; groups should separate along PC1/PC2 if the effect is strong.
- **Heatmap**: clear block structure by group indicates consistent DE patterns.
- **MA/Volcano**: many low-FDR points suggest robust signal; check symmetry and direction vs. expectations (e.g., P starvation up-regulates phosphate transport).


**Common red flags & quick fixes**  
- **Almost no DE genes** → Check replicate concordance (PCA), library sizes, or too-strict filters; consider reviewing strandedness/counting.
- **Weird direction** (e.g., expected transporter genes are “down”) → Make sure you’re reading `A_vs_B` correctly; “up” is **A over B**.
- **Gene IDs look unfamiliar in enrichment** → Ensure your enrichment mapping (`gene2go.tsv`) matches the same gene identifiers used in `DE_results`.

### 5) Functional enrichment of DEGs
#### What are we doing and why?
We’ll test whether DEGs are **over-represented** in **biological processes/pathways** (e.g., phosphate transport, root morphogenesis, ABA signaling). This connects statistics to **biology** and helps prioritize candidates that sit inside agronomically relevant pathways.

#### Option A — GOseq via Trinity 
Create with `vi`:
```bash
vi 05_trinity_goseq.sh
```
Press **`i`**, paste, then **`Esc` → `:wq`**:

```bash
#!/usr/bin/env bash
#SBATCH -p standard
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:25:00
#SBATCH -J day9_GOseq
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail

module purge
module load Trinity/2.15.2-foss-2023a
export TRINITY_HOME=${TRINITY_HOME:-$EBROOTTRINITY}
TRINITY_DE="$TRINITY_HOME/Analysis/DifferentialExpression"

# Choose your contrast name produced by Trinity (adjust as needed)
CONTR=Pminus_vs_Pplus

UP="edger_trinity/edgeR.DE_results/${CONTR}.DE_up.genes"
DOWN="edger_trinity/edgeR.DE_results/${CONTR}.DE_down.genes"
G2GO="metadata/gene2go.tsv"          # gene_id \t GO:xxxx;GO:yyyy
GLEN="counts/gene_lengths.txt"       # produced above from featureCounts

$TRINITY_DE/run_GOseq.pl \
  --genes_single_factor "$UP" \
  --gene2GO "$G2GO" \
  --lengths "$GLEN" \
  --out_prefix results_GOseq_up

$TRINITY_DE/run_GOseq.pl \
  --genes_single_factor "$DOWN" \
  --gene2GO "$G2GO" \
  --lengths "$GLEN" \
  --out_prefix results_GOseq_down
```
Submit:
```bash
sbatch 05_trinity_goseq.sh
```

### GOseq — Quick Evaluation (View Terms + Plot)
Minimal steps: **see enriched GO terms** (FDR ≤ 0.05) and **plot top terms**. Works for UP or DOWN outputs.

**1) View enriched GO terms (FDR ≤ 0.05)**  

**Pick the file** you want to inspect (replace with your actual filename):
```bash
# Examples (adjust to your real filenames produced by run_GOseq.pl)
UPTAB=results_GOseq_up.enriched.tsv
DOWNTAB=results_GOseq_down.enriched.tsv
```

**Print the top 20 enriched terms by FDR (header-aware):**  
```bash
TAB="$UPTAB"   # or: TAB="$DOWNTAB"

# Find the FDR column automatically and show the 20 best terms
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if(tolower($i) ~ /fdr|adj/) f=i; print; next}
            NR>1 && f>0 && $f!~"NA" {print}' "$TAB" \
| sort -t$'\t' -k${f},${f}g \
| head -n 20 \
| column -t
```

> If sorting fails because your shell doesn’t expand `${f}`, just open the file and skim:
```bash
column -t "$UPTAB" | less -S
```

**2) Plot top GO terms (simple barplot in R)**  

Create a tiny script **on the fly** and run it:
```bash
TAB="$UPTAB"   # or: TAB="$DOWNTAB"
Rscript - <<'RS'
tabfile <- Sys.getenv("TAB", "results_GOseq_up.enriched.tsv")
x <- read.table(tabfile, header=TRUE, sep="\t", quote="", comment.char="")
# Try to detect FDR and term/description columns
fdrcol <- grep("fdr|adj", tolower(colnames(x)), value=FALSE)[1]
termcol <- grep("term|description|name", tolower(colnames(x)), value=FALSE)[1]
if(is.na(fdrcol) || is.na(termcol)) stop("Could not find FDR or term column; open the file to check headers.")
x <- x[!is.na(x[,fdrcol]), ]
x$score <- -log10(pmax(x[,fdrcol], 1e-300))
o <- order(x$score, decreasing=TRUE)
x <- x[o, ]
topn <- head(x, 15)
pngfile <- sub("\\.\\w+$","",tabfile)
png(paste0(pngfile, "_top15.png"), width=900, height=600)
par(mar=c(5,20,2,2))
barplot(rev(topn$score), horiz=TRUE, names.arg=rev(topn[,termcol]), las=1,
        xlab="-log10(FDR)")
mtext("GOseq enrichment (top 15)", side=3, line=0.5)
dev.off()
cat("Wrote: ", paste0(pngfile, "_top15.png"), "\n")
RS
```

The script writes a PNG next to your table, e.g. `results_GOseq_up.enriched_top15.png`. Copy it to your laptop if needed:
```bash
scp USER@login.cluster:/work/USER/day09_rnaseq/results_GOseq_up.enriched_top15.png .
```

#### Option B — g:Profiler (Optional)
If species-specific GO is unavailable, we can demo with *Arabidopsis*:
```r
# Run interactively or via Rscript
library(gprofiler2)
up <- readLines('edger_trinity/edgeR.DE_results/Pminus_vs_Pplus.DE_up.genes')
res <- gost(query = up, organism = 'athaliana', correction_method = 'fdr')
write.csv(res$result, 'gprofiler_up.csv', row.names = FALSE)
```
> If you have an *Elaeis guineensis* gene2go, prefer **Option A** to stay species‑specific.

> *Breeding note:* Enrichment pinpoints processes (e.g., **phosphate transport**, **root development**, **ABA signaling**) to prioritize for selection and validation.

---
## RNAseq Troubleshooting
- **STAR mapping slow:** ensure `--readFilesCommand zcat` for gz FASTQs; check I/O throttling.
- **Few reads counted:** verify `-s` (strandedness) in featureCounts; wrong setting can halve counts.
- **No DE genes:** check replicate concordance (PCA), batch effects, and library sizes; consider increasing `--min_cpm` leniency for demo.
- **GOseq errors:** confirm `gene2go.tsv` format and `gene_lengths.txt` paths.

---

**Conclusion:**

This two-week intensive bioinformatics training course has provided a comprehensive overview of essential concepts and practical skills in handling and analyzing genomic data. From mastering the Linux command line to performing advanced RNA-seq analysis, you have gained valuable experience that will serve as a strong foundation for your future endeavors in bioinformatics. I encourage you to continue exploring, practicing, and applying these skills to unlock new biological insights. The field of bioinformatics is constantly evolving, and continuous learning is key to staying at the forefront of discovery. I wish you all the best in your bioinformatics journey!

# Congratulations!!! You have completed the training

---

### References
- Kong, L. et al. (2021). Comparative transcriptome analysis reveals novel insights into transcriptional responses to phosphorus starvation in oil palm root. *BMC Genomic Data*. https://pmc.ncbi.nlm.nih.gov/articles/PMC7863428/
- STAR — https://github.com/alexdobin/STAR
- Subread/featureCounts — http://subread.sourceforge.net/
- Trinity DE and GOseq wrappers — https://github.com/trinityrnaseq/trinityrnaseq/wiki
- edgeR — https://bioconductor.org/packages/edgeR
- g:Profiler — https://cran.r-project.org/package=gprofiler2

### Useful Tutorials and Resources

- [Introduction to RNAseq](https://scienceparkstudygroup.github.io/rna-seq-lesson/)  
- [RNAseq tutorial](https://sbc.shef.ac.uk/ngs_intro_workshop/03-rna-seq.nb.html)
- [Reference-based RNA-Seq data analysis](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html)
- [Informatics for RNA-seq](https://github.com/griffithlab/rnaseq_tutorial)
- [RNA-seq workflow](https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)

