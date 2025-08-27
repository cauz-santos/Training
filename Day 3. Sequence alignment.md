## Week 1: Foundations of Bioinformatics and Genomic Data Handling

### Day 3: Sequence Alignment and Dot Plot Analysis 

Day 3 will delve into the fundamental concepts of biological sequences and their comparison. Understanding DNA, RNA, and protein sequences is paramount in bioinformatics, as they are the raw data for most analyses. We will then explore how to visualize sequence similarities and differences using dot plot analysis, with a hands-on session in R.

**Morning Session: Basics of DNA, RNA, and Protein Sequences**

This session will provide a foundational understanding of the different types of biological sequences and their significance in molecular biology and genomics.

**DNA (Deoxyribonucleic Acid): The Blueprint of Life**

*   **Structure:** DNA is a double-stranded helix composed of nucleotides. Each nucleotide consists of a deoxyribose sugar, a phosphate group, and one of four nitrogenous bases: Adenine (A), Guanine (G), Cytosine (C), and Thymine (T).
*   **Base Pairing:** A always pairs with T, and C always pairs with G (Chargaff\'s rules). This complementary base pairing is crucial for DNA replication and repair.
*   **Function:** DNA carries the genetic instructions used in the growth, development, functioning, and reproduction of all known living organisms and many viruses. It serves as the long-term storage of genetic information.
*   **Representation in Bioinformatics:** DNA sequences are typically represented as strings of A, T, C, G characters.

**RNA (Ribonucleic Acid): The Messenger and More**

*   **Structure:** RNA is typically single-stranded and contains ribose sugar instead of deoxyribose. It has Uracil (U) instead of Thymine (T), so A pairs with U, and C pairs with G.
*   **Types and Functions:**
    *   **mRNA (messenger RNA):** Carries genetic information from DNA to ribosomes for protein synthesis.
    *   **tRNA (transfer RNA):** Carries specific amino acids to the ribosome during protein synthesis.
    *   **rRNA (ribosomal RNA):** A structural component of ribosomes, where protein synthesis occurs.
    *   **Non-coding RNAs (ncRNAs):** A diverse group of RNAs with regulatory and catalytic functions (e.g., microRNAs, long non-coding RNAs).
*   **Representation in Bioinformatics:** RNA sequences are represented as strings of A, U, C, G characters.

**Proteins: The Workhorses of the Cell**

*   **Structure:** Proteins are complex macromolecules made up of chains of amino acids linked by peptide bonds. There are 20 common amino acids, each with a unique side chain.
*   **Levels of Structure:** Proteins fold into specific three-dimensional structures (primary, secondary, tertiary, and quaternary) that determine their function.
*   **Function:** Proteins perform a vast array of functions within organisms, including catalyzing metabolic reactions (enzymes), DNA replication, responding to stimuli, and transporting molecules.
*   **Representation in Bioinformatics:** Protein sequences are represented as strings of single-letter amino acid codes (e.g., A for Alanine, L for Leucine).

**The Central Dogma of Molecular Biology:**

DNA -> RNA -> Protein. This fundamental concept describes the flow of genetic information within a biological system. Bioinformatics plays a critical role in analyzing and understanding each step of this process.

![Dogma](https://www.genome.gov/sites/default/files/media/images/2022-05/Central-dogma.jpg)
___
## Afternoon Session: Dot Plot Analysis for Sequence Comparisons

Dot plot analysis is a simple yet powerful graphical method for comparing two biological sequences (DNA, RNA, or protein) to identify regions of similarity, repeats, inversions, and deletions. It provides a visual representation of sequence alignment without performing a full, computationally intensive alignment.

**How Dot Plots Work:**

A dot plot is a two-dimensional matrix where one sequence is plotted along the x-axis and the other along the y-axis. A dot is placed at coordinates (i, j) if the character at position `i` in the x-axis sequence matches the character at position `j` in the y-axis sequence. Regions of similarity appear as diagonal lines.

*   **Perfect Match:** A single, continuous diagonal line.
*   **Repeats:** Multiple parallel diagonal lines.
*   **Inversions:** Diagonal lines perpendicular to the main diagonal.
*   **Insertions/Deletions:** Shifts or breaks in the diagonal lines.

**Hands-on Practice in R:**

R is a powerful statistical programming language widely used in bioinformatics for data analysis and visualization. We will use R to create dot plots and interpret their patterns.

**Required R Packages:**

We will likely use packages such as `seqinr` or `Biostrings` for sequence manipulation and base R plotting functions, or `ggplot2` for more advanced visualizations.

