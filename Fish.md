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

# Load modules
module load bwa
module load samtools

# Define paths
FASTA="kmers_50bp.fasta"
ORIG_REF="/lisc/scratch/course/pgbiow/data/genomes/Elaeis_guinensis_genomic.fna"
REF_DIR="$HOME/ref"
REF="$REF_DIR/Elaeis_guinensis_genomic.fna"
OUTDIR="results"
PREFIX="kmers_mapped"

# Prepare output folders
mkdir -p $OUTDIR $REF_DIR

# Copy and index reference if not already done
if [ ! -f "$REF.bwt" ]; then
    echo "Copying reference genome..."
    cp $ORIG_REF $REF

    echo "Indexing reference genome..."
    bwa index $REF
else
    echo "Reference genome and index already exist."
fi

# Step 1: Alignment
echo "Running BWA alignment..."
bwa aln -t 4 $REF $FASTA > $OUTDIR/$PREFIX.sai

echo "Generating SAM file..."
bwa samse $REF $OUTDIR/$PREFIX.sai $FASTA > $OUTDIR/$PREFIX.sam

# Step 2: Convert SAM to sorted BAM
echo "Converting SAM to sorted BAM..."
samtools view -bS $OUTDIR/$PREFIX.sam | samtools sort -o $OUTDIR/$PREFIX.sorted.bam
samtools index $OUTDIR/$PREFIX.sorted.bam

# Step 3: Filter uniquely mapped reads on chr2
echo "Filtering for unique mappings on chr2..."
samtools view -h $OUTDIR/$PREFIX.sorted.bam chr2 | \
  awk '$0 ~ /^@/ || $5 >= 60' > $OUTDIR/${PREFIX}_chr2_unique.sam

# Step 4: Extract unique read IDs
echo "Extracting read IDs..."
awk '$1 !~ /^@/ {print $1}' $OUTDIR/${PREFIX}_chr2_unique.sam | sort | uniq > $OUTDIR/chr2_unique_ids.txt

# Step 5: Extract FASTA sequences from original input
echo "Extracting FASTA entries..."
awk 'NR==FNR {ids[$1]; next} /^>/ {header=$0; id=substr($0, 2); show=(id in ids)} show' \
    $OUTDIR/chr2_unique_ids.txt $FASTA > $OUTDIR/chr2_unique_kmers.fasta

echo "✅ Mapping and extraction complete. Results in $OUTDIR"
```
Submit to SLURM
sbatch run_bwa_map.sh
