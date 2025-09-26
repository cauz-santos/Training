Step 1: Convert probes.txt to FASTA

Script:  vi convert_probes_to_fasta.py

Save this script if you haven't already:

```bash
#convert_probes_to_fasta.py
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as infile, open(output_file, "w") as out:
    for line in infile:
        if line.strip() == "":
            continue
        parts = line.strip().split("\t")
        if len(parts) != 2:
            continue
        probe_id, seq = parts
        out.write(f">{probe_id}\n{seq.strip()}\n")
```

Run it manually on the cluster
```bash
python3 convert_probes_to_fasta.py probes.txt kmers_50bp.fasta
```


This will create kmers_50bp.fasta in your working directory.

Step 2: Run BWA using the FASTA file
Script: vi run_bwa_map.sh

```bash
#!/bin/bash
#SBATCH --job-name=bwa_map_kmers
#SBATCH --output=logs/bwa_map_%j.out
#SBATCH --error=logs/bwa_map_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G

# Load required modules
module load bwa
module load samtools
module load seqkit  # optional

# INPUTS
FASTA="kmers_50bp.fasta"                     # your converted FASTA
REF="/lisc/scratch/course/pgbiow/data/genomes/Elaeis_guinensis_genomic.fna"          # CHANGE this to your actual BWA index path
OUTDIR="results"
PREFIX="kmers_mapped"

mkdir -p $OUTDIR

echo "Aligning kmers..."
bwa aln -t 4 $REF $FASTA > $OUTDIR/$PREFIX.sai

echo "Generating SAM..."
bwa samse $REF $OUTDIR/$PREFIX.sai $FASTA > $OUTDIR/$PREFIX.sam

echo "Sorting and indexing..."
samtools view -bS $OUTDIR/$PREFIX.sam | samtools sort -o $OUTDIR/$PREFIX.sorted.bam
samtools index $OUTDIR/$PREFIX.sorted.bam

echo "Filtering uniquely mapped kmers on chr2 (MAPQ ≥ 60)..."
samtools view -h $OUTDIR/$PREFIX.sorted.bam chr2 | \
  awk '$0 ~ /^@/ || $5 >= 60' > $OUTDIR/${PREFIX}_chr2_unique.sam

echo "Extracting unique read IDs..."
awk '$1 !~ /^@/ {print $1}' $OUTDIR/${PREFIX}_chr2_unique.sam | sort | uniq > $OUTDIR/chr2_unique_ids.txt

if command -v seqkit &> /dev/null; then
    echo "Extracting FASTA entries for unique chr2 probes..."
    seqkit grep -f $OUTDIR/chr2_unique_ids.txt $FASTA > $OUTDIR/chr2_unique_kmers.fasta
else
    echo "seqkit not found – skipping FASTA extraction."
fi

echo  "Finished."
```
Submit to SLURM
sbatch run_bwa_map.sh
