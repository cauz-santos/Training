## Week 2: Applied Bioinformatics for Genomics and Breeding

### Day 6: SNP Calling, Variant Processing, and Phylogenetics

Week 2 shifts our focus to applied bioinformatics, specifically in the context of genomics and breeding. Day 6 will introduce you to Single Nucleotide Polymorphisms (SNPs), their detection (SNP calling), processing variant data, and how to use genetic variation to understand evolutionary relationships through phylogenetics.

**Morning Session 1: Understanding SNPs and DArT markers (RAN to provide an introduction of Verdant Tag Panel)**

**Single Nucleotide Polymorphisms (SNPs):**

SNPs are the most common type of genetic variation among individuals. A SNP occurs when a single nucleotide (A, T, C, or G) in the genome differs between members of a species or paired chromosomes in an individual. For example, at a specific position in the genome, one individual might have an 'A' while another has a 'G'.

*   **Significance:** SNPs are crucial genetic markers used in various applications:
    *   **Disease Association:** Identifying SNPs associated with disease susceptibility or resistance.
    *   **Pharmacogenomics:** Predicting an individual's response to drugs.
    *   **Population Genetics:** Studying genetic diversity and evolutionary history.
    *   **Plant and Animal Breeding:** Marker-assisted selection for desirable traits.
*   **Frequency:** To be considered a SNP, the variation must be present in at least 1% of the population.

**DArT (Diversity Arrays Technology) markers:**

DArT is a high-throughput genotyping platform used to discover and score thousands of genetic markers (often SNPs and presence/absence variations) across a genome. It is particularly useful for species with complex genomes or without a reference genome. DArT markers are often used in breeding programs for genetic mapping, diversity analysis, and marker-assisted selection.

*   **Verdant Tag Panel (RAN to provide introduction):** This refers to a specific DArT marker panel, likely optimized for certain plant species or breeding objectives. The trainer (RAN) will provide specific details on its application and utility.

**Morning Session 2: SNP Calling and Annotation (e.g., GATK, Stacks, bcftools, VCFtools, SnpEff)**

SNP calling is the process of identifying genetic variations (SNPs and small insertions/deletions, or indels) from sequencing data. This typically involves aligning sequencing reads to a reference genome and then identifying positions where the observed nucleotides differ from the reference.

**General SNP Calling Workflow:**

1.  **Read Alignment:** Align raw sequencing reads (FASTQ) to a reference genome (e.g., using BWA, Bowtie2) to produce BAM files.
2.  **Variant Discovery:** Identify potential variant sites based on discrepancies between reads and the reference.
3.  **Variant Genotyping:** Determine the genotype (e.g., homozygous reference, heterozygous, homozygous alternative) at each variant site for each sample.
4.  **Variant Filtering:** Apply filters to remove low-quality or artifactual variants.

**Common Tools for SNP Calling:**

    *   **GATK (Genome Analysis Toolkit):** A widely used, industry-standard suite of tools for variant discovery and genotyping, particularly for human genome data. It includes best practices workflows for germline and somatic variant calling.
        ```bash
        # Example GATK workflow (conceptual, requires Java and GATK installed)
        # 1. Index reference genome (if not already done)
        # samtools faidx reference.fasta
        # java -jar picard.jar CreateSequenceDictionary R=reference.fasta O=reference.dict

        # 2. HaplotypeCaller for variant discovery (per sample)
        # gatk HaplotypeCaller -R reference.fasta -I aligned_reads.bam -O raw_variants.g.vcf.gz -ERC GVCF

        # 3. CombineGVCFs (if multiple samples)
        # gatk CombineGVCFs -R reference.fasta -V sample1.g.vcf.gz -V sample2.g.vcf.gz -O combined.g.vcf.gz

        # 4. GenotypeGVCFs (joint genotyping)
        # gatk GenotypeGVCFs -R reference.fasta -V combined.g.vcf.gz -O raw_variants.vcf.gz

        # 5. Variant Quality Score Recalibration (VQSR) - for large datasets
        # gatk VariantRecalibrator -R reference.fasta -V raw_variants.vcf.gz ... -O recal.tranches
        # gatk ApplyVQSR -R reference.fasta -V raw_variants.vcf.gz ... -O final_variants.vcf.gz
        ```
    *   **Stacks:** A software pipeline for building loci from short-read sequences, particularly useful for SNP discovery and genotyping in populations without a reference genome (de novo assembly) or with a reference (reference-based assembly).
        ```bash
        # Example Stacks workflow (conceptual)
        # 1. Process_radtags (demultiplexing and quality filtering)
        # process_radtags -f raw_reads.fastq -o ./output_dir -b barcodes.txt -e sbfI -r -c -q -t 90 --inline_null

        # 2. Ustacks (de novo assembly of loci within individuals)
        # ustacks -f sample_1.fq -o ./output_dir -i 1 -m 3 -M 2 -p 4

        # 3. Cstacks (build catalog of loci across individuals)
        # cstacks -b 1 -P ./output_dir -s sample1 -s sample2 -s sampleN

        # 4. Sstacks (match samples to catalog)
        # sstacks -b 1 -P ./output_dir -s sample1 -s sample2 -s sampleN

        # 5. Gstacks (genotype individuals from aligned reads)
        # gstacks -b 1 -P ./output_dir -t 4

        # 6. Populations (filter and export VCF)
        # populations -P ./output_dir -M popmap.txt -r 0.8 --vcf
        ```
    *   **bcftools:** A suite of command-line utilities for manipulating and analyzing VCF (Variant Call Format) files. It can be used for filtering, merging, and comparing VCF files, as well as for basic variant calling.
        ```bash
        # Example bcftools variant calling workflow
        # 1. Align reads (e.g., with BWA) and sort/index BAM (as shown in Hands-on section)

        # 2. Call variants using mpileup and call
        bcftools mpileup -Ou -f reference.fasta aligned_reads.sorted.bam | bcftools call -mv -Ov -o variants.vcf

        # 3. Filter variants (example: remove low quality variants)
        bcftools view -i 'QUAL>20 && FORMAT/DP>10' variants.vcf -o filtered_variants.vcf
        ```
    *   **VCFtools:** Another set of command-line tools for working with VCF files, offering functionalities for filtering, statistics, and format conversion.
        ```bash
        # Example VCFtools usage (filtering and statistics)
        # Filter variants by minor allele frequency (MAF) and missing data
        vcftools --vcf variants.vcf --maf 0.05 --max-missing 0.8 --recode --recode-INFO-all --out filtered_by_maf_missing

        # Calculate SNP diversity (pi)
        vcftools --vcf filtered_variants.vcf --site-pi --out site_pi
        ```

**SNP Annotation:**

Once SNPs are called, they need to be annotated to understand their potential functional impact. Annotation involves determining if a SNP is in a gene, what type of change it causes (e.g., synonymous, non-synonymous, stop-gain), and if it has been previously reported in databases.

*   **SnpEff:** A popular tool for annotating genetic variants. It predicts the effect of variants on genes and proteins (e.g., missense, nonsense, frameshift) and can integrate with various databases.
    ```bash
    # Example SnpEff workflow
    # 1. Download and build the SnpEff database for your organism (e.g., for human GRCh38)
    # java -jar snpEff.jar download GRCh38

    # For custom genome/annotation (requires GTF/GFF3 and FASTA)
    # Create a snpEff.config file pointing to your genome and annotation
    # java -jar snpEff.jar build -gff3 -v your_organism_database

    # 2. Annotate your VCF file
    java -jar snpEff.jar your_organism_database filtered_variants.vcf > annotated_variants.vcf
    ```

**Afternoon Session 3: Hands-on Practice: SNP Calling Workflow**

In this hands-on session, we will perform a simplified SNP calling workflow using a small dataset. We will focus on practical application of some of the tools discussed.

**Example Workflow Steps (using a subset of tools):**

1.  **Reference Genome Preparation:** Indexing the reference genome for alignment.
    ```bash
    bwa index reference.fasta
    samtools faidx reference.fasta
    java -jar picard.jar CreateSequenceDictionary R=reference.fasta O=reference.dict
    ```
2.  **Read Alignment (e.g., BWA-MEM):**
    ```bash
    bwa mem reference.fasta raw_reads_R1.fastq.gz raw_reads_R2.fastq.gz | samtools view -bS - > aligned_reads.bam
    ```
3.  **Sort and Index BAM:**
    ```bash
    samtools sort aligned_reads.bam -o aligned_reads.sorted.bam
    samtools index aligned_reads.sorted.bam
    ```
4.  **Variant Calling (e.g., bcftools):**
    ```bash
    bcftools mpileup -Ou -f reference.fasta aligned_reads.sorted.bam | bcftools call -mv -Ov -o variants.vcf
    ```
5.  **Variant Filtering (e.g., VCFtools):**
    ```bash
    vcftools --vcf variants.vcf --remove-filtered-all --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out filtered_variants
    ```
6.  **Variant Annotation (e.g., SnpEff):**
    ```bash
    # First, download and build the SnpEff database for your organism
    # java -jar snpEff.jar download GRCh38 # Example for human
    java -jar snpEff.jar build -gff3 -v your_organism_database # For custom database
    java -jar snpEff.jar your_organism_database filtered_variants.recode.vcf > annotated_variants.vcf
    ```

**Afternoon Session 4: Hands-on: Constructing Phylogenetic Trees**

Phylogenetics is the study of evolutionary relationships among groups of organisms (or genes) through time. Phylogenetic trees are diagrams that depict these relationships. They are constructed based on genetic (or morphological) data, often using sequence variations like SNPs.

**Key Concepts in Phylogenetics:**

*   **Taxa/OTUs:** The operational taxonomic units (e.g., species, individuals, genes) at the tips of the branches.
*   **Nodes:** Represent common ancestors.
*   **Branches:** Represent evolutionary lineages.
*   **Root:** The common ancestor of all taxa in the tree (if rooted).
*   **Clade:** A group of organisms that includes an ancestor and all of its descendants.

**Methods for Tree Construction:**

*   **Distance-based methods (e.g., Neighbor-Joining):** Calculate genetic distances between sequences and then build a tree based on these distances.
*   **Parsimony methods:** Find the tree that requires the fewest evolutionary changes (mutations) to explain the observed data.
*   **Maximum Likelihood (ML) methods:** Find the tree that maximizes the probability of observing the given data under a specific evolutionary model.
*   **Bayesian methods:** Use Bayesian inference to estimate the posterior probability of trees.

**Hands-on Practice with Phylogenetic Tools:**

We will use a set of aligned sequences (e.g., from a multi-sample VCF converted to FASTA) and apply tools to construct and visualize phylogenetic trees.

**Tools (examples):**

    *   **MEGA (Molecular Evolutionary Genetics Analysis):** A user-friendly software for sequence alignment, phylogenetic tree construction, and evolutionary analysis. (Often GUI-based, but command-line versions or specific modules can be used).
    *   **RAxML/IQ-TREE:** Popular command-line tools for maximum likelihood phylogenetic inference, known for their speed and accuracy.
        ```bash
        # Example RAxML-NG workflow (conceptual)
        # 1. Prepare input alignment (e.g., FASTA format)
        # 2. Run RAxML-NG for tree inference and bootstrapping
        raxml-ng --msa input_alignment.fasta --model GTR+G --prefix T1 --threads 4 --bs-trees 100
        # --msa: input multiple sequence alignment
        # --model: evolutionary model (e.g., GTR+G for nucleotides)
        # --prefix: output file prefix
        # --threads: number of CPU threads
        # --bs-trees: number of bootstrap replicates
        ```
    *   **FastTree:** A fast algorithm for constructing approximate maximum-likelihood phylogenetic trees from large alignments.
        ```bash
        # Example with FastTree (assuming input.fasta is your aligned sequence file)
        FastTree -gtr -nt input.fasta > tree.nwk
        # -gtr: General Time Reversible model (for nucleotide data)
        # -nt: Specifies nucleotide data
        # tree.nwk: Output tree in Newick format
        ```
    *   **FigTree/iTOL:** Tools for visualizing and annotating phylogenetic trees.

**Example Workflow (Conceptual):**

1.  **Prepare Input Data:** Convert your VCF file (containing SNPs from multiple samples) into a format suitable for phylogenetic analysis (e.g., FASTA alignment). This often involves custom scripting or tools like `vcf2phylip`.
    ```bash
    # Conceptual step: Convert VCF to FASTA alignment for phylogenetic analysis
    # This often requires a custom script or tool like vcf2phylip.py
    # python vcf2phylip.py -i your_variants.vcf -o alignment.phy
    # Or for simple FASTA from VCF (for SNP sites only):
    # bcftools query -f ">[%SAMPLE\t%REF%ALT]\n" your_variants.vcf | sed 's/\t//g' > alignment.fasta
    ```
2.  **Run Tree Construction:** Execute a phylogenetic program.
    ```bash
    # Using FastTree for a quick tree
    FastTree -gtr -nt alignment.fasta > tree.nwk
    ```
3.  **Visualize Tree:** Open the `.nwk` file in a tree visualization software like FigTree or upload to iTOL.

By the end of Day 6, you will have a practical understanding of how to identify genetic variations and use them to infer evolutionary relationships, skills critical for both basic research and applied breeding programs.


