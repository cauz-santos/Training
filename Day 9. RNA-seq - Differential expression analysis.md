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

| Group | Time | Tissue | Rep | SRA Run | Suggested sample name |
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


### 1) STAR genome index (Slurm)

**What are we doing and why?**  
We will **build a STAR genome index** from a reference **FASTA** and **GTF**. STAR needs this index to align reads quickly and to detect **splice junctions** accurately.  
**Breeding relevance:** robust splice-aware alignment preserves signals from genes and isoforms involved in **phosphorus-use efficiency (PUE)** and other stress responses—exactly the biology we want to prioritize for selection, marker design, and follow-up validation.

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
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 00:25:00
#SBATCH -J day9_star_index
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail

module purge
module load Trinity/2.15.2-foss-2023a
module load STAR/2.7.10a

REF=ref/reference.fa
GTF=ref/annotation.gtf
IDX=ref/STAR_index
READLEN=${READLEN:-100}          # set to your read length (100 or 150)
SJDB=$((READLEN-1))

mkdir -p "$IDX"
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir "$IDX" \
     --genomeFastaFiles "$REF" \
     --sjdbGTFfile "$GTF" \
     --sjdbOverhang $SJDB
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
**Breeding relevance:** consistent, high-quality alignment across treatments (−P vs +P) is essential to avoid false differential expression and to capture splice-aware signals in genes underlying **phosphorus-use efficiency (PUE)**.

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

IDX=ref/STAR_index
SAMPLE=${SAMPLE:-Pminus_28d_rep1}
R1=fastq/${SAMPLE}_R1.fastq.gz
R2=fastq/${SAMPLE}_R2.fastq.gz
OUT=star/${SAMPLE}
mkdir -p "$OUT"

STAR --runThreadN 16 \
  --genomeDir "$IDX" \
  --readFilesIn "$R1" "$R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUT/" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS nM XS \
  --twopassMode Basic

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
**Breeding relevance:** accurate gene counts are the backbone of **edgeR** differential expression. This is where nutrient/stress pathways (e.g., phosphate transport, root development) first become **quantifiable** for ranking candidates and informing selection decisions.

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
module load Subread/2.0.3

GTF=ref/annotation.gtf
OUT=counts
mkdir -p "$OUT"

# collect all BAMs (pre-mapped + demo)
ls star/*/Aligned.sortedByCoord.out.bam > ${OUT}/bam_list.txt

# set strandedness: 0=unstranded, 1=forward, 2=reverse
STRAND=${STRAND:-0}

featureCounts -T 16 -p -B -C -s ${STRAND} \
  -t exon -g gene_id \
  -a "$GTF" \
  -o ${OUT}/featureCounts.txt \
  $(cat ${OUT}/bam_list.txt)

# Clean genes x samples matrix
python - << 'PY'
import pandas as pd
fc = pd.read_csv('counts/featureCounts.txt', sep='\\t', comment='#')
cts = fc.iloc[:, [0] + list(range(6, fc.shape[1]))].copy()
cts.rename(columns={'Geneid':'gene_id'}, inplace=True)
cts.to_csv('counts/counts_matrix.tsv', sep='\\t', index=False)
# gene lengths for GOseq
fc[['Geneid','Length']].rename(columns={'Geneid':'gene_id','Length':'length'}) \
  .to_csv('counts/gene_lengths.txt', sep='\\t', index=False, header=False)
PY
```

Submit:
```bash
# Example: reverse-stranded libraries
STRAND=2 sbatch 03_featurecounts.sh
```

> *Breeding note:* Gene‑level counts quantify **stress and nutrient pathways** driving resilience and yield stability under low‑P.

### After featureCounts: how to inspect your outputs

#### 1) Inspect the summary table
The file `featureCounts.txt.summary` reports, **per sample**, how many fragments were `Assigned` (counted to a gene) and how many were `Unassigned_*` for various reasons.

```bash
# Quick view
less -S featureCounts.txt.summary

# Nicely aligned columns (press q to quit)
column -t featureCounts.txt.summary | less -S
```
**What to look for:**
- `Assigned` should be the **largest** row across samples.
- Large `Unassigned_NoFeatures` → GTF mismatch (wrong genome/annotation).
- Large `Unassigned_Strand` → wrong `-s` value (strandedness).
- Large `Unassigned_MultiMapping` → repetitive regions or liberal aligner settings.

**Compute Assigned % per sample (quick check):**
```bash
awk 'NR==1{for(i=2;i<=NF;i++) S[i]=$i; next} {for(i=2;i<=NF;i++) T[i]+=$i} $1=="Assigned"{for(i=2;i<=NF;i++) A[i]=$i}
     END{for(i=2;i<=NF;i++){pct=(A[i]/T[i])*100; printf "%s\t%.2f%% Assigned\n", S[i], pct}}' featureCounts.txt.summary
```

#### 2) Spot-check the clean counts matrix
The file `counts/counts_matrix.tsv` is **genes × samples** and is what you’ll feed into Trinity/edgeR.

```bash
# Header line (gene_id and sample names, check order!)
head -n 1 counts_matrix.tsv

# First 10 data rows
head counts_matrix.tsv | column -t
```

**Check total library sizes (sum of counts) per sample:**
```bash
awk -F'\t' 'NR==1{for(i=2;i<=NF;i++) name[i]=$i; next}
            {for(i=2;i<=NF;i++) lib[i]+=$i}
            END{for(i=2;i<=NF;i++) printf "%s\t%d\n", name[i], lib[i]}' counts_matrix.tsv | column -t
```

**Sanity check for all-zero genes (should be few or none):**
```bash
awk -F'\t' 'NR>1{sum=0; for(i=2;i<=NF;i++) sum+=$i; if(sum==0) print $1}' counts_matrix.tsv | head
```

**Ensure sample order matches Trinity’s `samples.txt`:**
```bash
echo "Counts header:"
head -n1 counts_matrix.tsv

echo -e "\nTrinity groups/samples:"
cat metadata/samples.txt
```

**Common red flags & fixes**
- **Assigned < 60%** across samples → check `-s` and GTF/genome match.
- **One sample much smaller library size** → expect more dispersion; check `Log.final.out` and earlier QC.
- **Header/sample name mismatch** → rename BAMs consistently or regenerate the matrix with correct order.


### 4) Differential expression with **edgeR via Trinity** (Slurm)
Create `metadata/samples.txt` as shown above. Then create `scripts/04_trinity_edger.sh`:
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

module purge
module load Trinity/2.15.2-foss-2023a
export TRINITY_HOME=${TRINITY_HOME:-$EBROOTTRINITY}
TRINITY_DE="$TRINITY_HOME/Analysis/DifferentialExpression"

OUTDE=edger_trinity
mkdir -p "$OUTDE"

$TRINITY_DE/run_DE_analysis.pl \
  --matrix counts/counts_matrix.tsv \
  --method edgeR \
  --samples_file metadata/samples.txt \
  --output "$OUTDE" \
  --min_cpm 1 --min_reps_min_cpm 2

# Optional: clustering/PCA/heatmaps
$TRINITY_DE/analyze_diff_expr.pl \
  --matrix ${OUTDE}/counts_TMM_normalized.matrix \
  --samples metadata/samples.txt \
  --min_rowSum_counts 10 \
  --P 0.05 --C 1 \
  --output ${OUTDE}/diffexpr_plots
```
Submit:
```bash
sbatch scripts/04_trinity_edger.sh
```

**Outputs of interest**
- `edger_trinity/edgeR.DE_results/*DE_results` (per-contrast gene tables with logFC, FDR)
- `edger_trinity/diffexpr_plots/` (PCA, heatmaps, volcano/MA)

> *Breeding note:* DE highlights candidate genes/TFs for **marker design**, **introgression**, or **editing** to improve PUE.

---

## 5) Functional enrichment of DEGs
### Option A — GOseq via Trinity (length bias aware)
Create `scripts/05_trinity_goseq.sh`:
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
sbatch scripts/05_trinity_goseq.sh
```

### Option B — g:Profiler (quick)
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

## Troubleshooting
- **STAR mapping slow:** ensure `--readFilesCommand zcat` for gz FASTQs; check I/O throttling.
- **Few reads counted:** verify `-s` (strandedness) in featureCounts; wrong setting can halve counts.
- **No DE genes:** check replicate concordance (PCA), batch effects, and library sizes; consider increasing `--min_cpm` leniency for demo.
- **GOseq errors:** confirm `gene2go.tsv` format and `gene_lengths.txt` paths.

---

**Conclusion:**

This two-week intensive bioinformatics training course has provided a comprehensive overview of essential concepts and practical skills in handling and analyzing genomic data. From mastering the Linux command line to performing advanced RNA-seq analysis, you have gained valuable experience that will serve as a strong foundation for your future endeavors in bioinformatics. I encourage you to continue exploring, practicing, and applying these skills to unlock new biological insights. The field of bioinformatics is constantly evolving, and continuous learning is key to staying at the forefront of discovery. I wish you all the best in your bioinformatics journey!

# Congratulations!!! You have completed the training

---

## References
- Kong, L. et al. (2021). Comparative transcriptome analysis reveals novel insights into transcriptional responses to phosphorus starvation in oil palm root. *BMC Genomic Data*. https://pmc.ncbi.nlm.nih.gov/articles/PMC7863428/
- STAR — https://github.com/alexdobin/STAR
- Subread/featureCounts — http://subread.sourceforge.net/
- Trinity DE and GOseq wrappers — https://github.com/trinityrnaseq/trinityrnaseq/wiki
- edgeR — https://bioconductor.org/packages/edgeR
- g:Profiler — https://cran.r-project.org/package=gprofiler2
