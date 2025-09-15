## Week 1: Foundations of Bioinformatics and Genomic Data Handling

## Day 2: Bioinformatics Scripting Essentials (Bash, Python)

On Day 2, we transition from using individual commands to combining them into scripts. Scripting is the key to automating repetitive tasks, creating reproducible workflows, and handling large datasets efficiently. We will cover the fundamentals of scripting in two powerful languages for bioinformatics: Bash and Python.

## Introduction to Bash Scripting

Bash (Bourne-Again SHell) is the default command-line interpreter on most Linux systems. A Bash script is simply a text file containing a series of commands that are executed sequentially. It's an excellent tool for automating command-line workflows.

**Creating a Simple Bash Script:**

1.  **Create a script file:** Use a text editor like `nano` or `vim` to create a file, for example, `my_first_script.sh`.
2.  **Add the shebang:** The first line of a Bash script should be `#!/bin/bash`. This tells the system to use the Bash interpreter to run the script.
3.  **Add commands:** Add the commands you want to execute, one per line.
4.  **Make the script executable:** Use the `chmod +x my_first_script.sh` command to give the script execute permissions.
5.  **Run the script:** Execute the script by typing `./my_first_script.sh`.


### 1) Example of Basic Bash Script:

**Step 1: Create a script file**   

We first need to create a new file that will hold our script.    
Use the `nano` text editor (you can also use `vim` or any other editor):

```bash
nano hello.sh
```

**Step 2: Write the script**  

Inside the editor, type the following lines:
```bash
#!/bin/bash
# This is my first Bash script
# Lines starting with # are comments and not executed

echo "Hello, Bioinformatics World!"
```

When finished, save and exit:

In `nano`: press `Ctrl + O` → `Enter` → `Ctrl + X`. 

**Step 3: Make the script executable**  

By default, new files are not executable. We need to give it execute permission:  
```bash
chmod +x hello.sh
```

**Step 4: Run the script**  

Now we can run it by calling the file directly:  
```bash
./hello.sh
```

### 2) Example Bash Script for Automation:

**Step 1 – Create directories**
Let's create a script that automates the process of creating a new project directory with subdirectories for data, scripts, and results.

Please create a file typing `vi create_project.sh`in the terminal, and then copy the following content:

```bash
#!/bin/bash

# ============================================
# Bioinformatics Project Structure Creator
# ============================================
# Usage:
#   ./create_project.sh <project_name>
#
# Example:
#   ./create_project.sh my_bio_project
# ============================================

# Check if the user provided a name
if [ -z "$1" ]; then
  echo "Error: No project name supplied."
  echo "Usage: ./create_project.sh <project_name>"
  exit 1
fi

PROJECT_NAME=$1

echo "Creating project: $PROJECT_NAME"

mkdir -p "$PROJECT_NAME"/{raw_data,processed_data,scripts,results,logs}

echo "Structure created:"
tree "$PROJECT_NAME"
```

To run this script, you would save it as `create_project.sh`, make it executable (`chmod +x create_project.sh`), and then run it with a project name as an argument: `./create_project.sh my_new_project`.

For example:
```bash
./create_project.sh Test
```

**Step 2 – Rename Files from _raw.txt to _cleaned.txt**
Let’s rename all raw input files to mark them as cleaned. Please create a file typing `vi rename_files.sh`in the terminal, and then copy and paste the following content:

```bash
#!/bin/bash

for file in *_raw.txt; do
    newname=$(echo "$file" | sed 's/_raw/_cleaned/')
    mv "$file" "$newname"
    echo "Renamed $file → $newname"
done
```

To run this script, you would save it and make it executable (`chmod +x rename_files.sh`), and then run it typing in the terminal: `./rename_files.sh`.

**Step 3 – Count Files Containing “Status: OK” or “Status: FAIL”**
Let’s analyze the content of the cleaned files and summarize how many contain "Status: OK" or "Status: FAIL"  
Create the file script using `vi summarize_status.sh`

```bash
#!/bin/bash

# ============================================
# Status Summary Report
# ============================================

OK_COUNT=0
FAIL_COUNT=0

for file in *_cleaned.txt; do
    if grep -q "Status: OK" "$file"; then
        OK_COUNT=$((OK_COUNT + 1))
    elif grep -q "Status: FAIL" "$file"; then
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
done

echo "Summary Report:"
echo "Files with Status: OK → $OK_COUNT"
echo "Files with Status: FAIL → $FAIL_COUNT"
```

To run this script, you would save it and make it executable (`chmod +x summarize_status.sh`), and then run it typing in the terminal: `./summarize_status.sh`.

### Exercise: Write a Bash Script to Summarize a Genome FASTA

**Goal:** Learn to write a Bash script that processes the *Elaeis guineensis* reference genome FASTA and reports some simple statistics.


### What you will do
- Create a new Bash script called `genome_stats.sh`.
- Make it executable and run it on the provided genome FASTA.
- Learn how to **combine commands** into a reusable workflow.


**Step 1: Create the script file**  

Open a new script file:

```bash
vi genome_stats.sh
```

**Step 2: Add this content**  
```bash
#!/bin/bash
# Simple Bash script to summarize a genome FASTA file

# Usage: ./genome_stats.sh <genome.fasta>

FASTA=$1

echo "Analyzing genome file: $FASTA"
echo "-----------------------------------"

# Count the number of sequences (headers start with ">")
echo -n "Number of sequences: "
grep -c "^>" "$FASTA"

# Count the total number of bases (ignore header lines)
echo -n "Total bases: "
grep -v "^>" "$FASTA" | tr -d '\n' | wc -c

# Report the first 3 sequence names
echo "First 3 sequence IDs:"
grep "^>" "$FASTA" | head -n 3
```

**Step 3: Make it executable**
```bash
chmod +x genome_stats.sh
```

**Step 4: Run the script**
```bash
./genome_stats.sh data/Elaeis_guineensis.fasta
```
*Please check the basic statistics from the output  

___
## Extra: Useful Unix Tools (Pipes, head/tail, grep, cut)

Before moving to Python, let’s explore a few essential Unix tools that make working with text files powerful and efficient.  
These are often used inside **pipes (`|`)**, which pass the output of one command as input to the next.


**1. Pipes ( | )**  
Pipes allow you to **combine commands** into a single workflow.

```bash
# Count how many .txt files are in the folder
ls *.txt | wc -l
```
- ls *.txt lists all .txt files.
- wc -l counts the number of lines.

**2. head and tail**  
`head` shows the first lines of a file, while `tail` shows the last lines.

```bash
# Show the first 5 lines of a table
head -5 sample_metadata.tsv

# Show the last 5 lines of the same file
tail -5 sample_metadata.tsv
```

**3. grep**  
`grep` searches for patterns in text files.

```bash
# Find all lines containing "Status: OK"
grep "Status: OK" *_cleaned.txt

# Count how many times "Status: FAIL" appears
grep -c "Status: FAIL" *_cleaned.txt
```
- grep "pattern" file prints lines containing that pattern.
- -c counts matches instead of printing them.

**4. cut**  
`cut` extracts columns from tabular data.  
By default, it assumes the file is tab-delimited (\t).

```bash
# Extract the second column (Treatment) from a metadata table
cut -f2 sample_metadata.tsv | head
```


___

## Introduction to Python for Bioinformatics

Python is a versatile, high-level programming language that has become a cornerstone of bioinformatics due to its readability, extensive libraries (like Biopython), and powerful data analysis capabilities. We will cover the absolute basics to get you started.

This section includes **two beginner-friendly exercises** to practice essential skills like:

- Lists and loops
- File reading and filtering
- Conditional logic
- String manipulation


**Table 1. Python Vocabulary: Key Concepts to Get Started**

| Concept             | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| `print()`           | Displays output on the screen.                                              |
| **Variable**        | A name used to store a value (e.g., `x = 5`).                               |
| **String**          | Text data, written in quotes (e.g., `"ATGCGT"`).                            |
| **Integer / Float** | Whole numbers (`int`) or decimals (`float`).                               |
| **List**            | An ordered collection of values (e.g., `[1, 2, 3]`).                        |
| **Dictionary**      | A collection of key-value pairs (e.g., `{"gene": "ATG"}`).                 |
| `for` loop          | Repeats an action for each item in a list or range.                        |
| `if` statement      | Executes code only if a condition is true.                                 |
| `def`               | Used to define a **function** — a reusable block of code.                  |
| `import`            | Loads external libraries or modules.                                       |
| `open()`            | Opens a file for reading or writing.                                       |
| `readline()` / `readlines()` | Reads data from a file.                                     |
| `len()`             | Returns the number of elements in a list or characters in a string.        |
| `range()`           | Generates a sequence of numbers, often used in loops.                      |
| `type()`            | Returns the data type of a variable.                                       |


### Open the Python Interpreter

In your terminal, type:

```bash
python3
```

You should see something like:
```bash
Python 3.10.12 (default, ...)
>>>
```
This >>> prompt means Python is ready for commands.


Let's type some basic commands:

Exercise 1: Print a Welcome Message  
Type the following commands in the Python prompt:

```bash
print("Welcome to Bioinformatics with Python!")
```

Exercise 2: Do Basic Math with Variables:
```bash
x = 5
y = 10
print("Sum:", x + y)
print("Product:", x * y)
```

Exercise 3: Create and Print a List:
```bash
samples = ["Sample1", "Sample2", "Sample3"]
print("Total samples:", len(samples))
print("First sample:", samples[0])
```

Exercise 5: Define and Use a Function:
```bash
def gc_content(sequence):
    gc = sequence.count("G") + sequence.count("C")
    return gc / len(sequence) * 100

seq = "ATGCGCGTAGCTAGC"
print("GC content:", gc_content(seq), "%")
```

Next, we explore some Python scripts for automating tasks...

### Example 1: Filtering a List of Gene Names  

**Goal:** Print only gene names that start with `"ABC"` from a predefined list.

**Concepts Covered:**
- Lists
- For loops
- String methods (`startswith`)
- Conditional statements

Open a file with `vi gene_filter.py`, copy and paste:

```python
# gene_filter.py

genes = ["ATG5", "ABC1", "NBS1", "ATP6", "COX3", "ABC2", "PPR1"]

print("Genes starting with 'ABC':")
for gene in genes:
    if gene.startswith("ABC"):
        print(f" - {gene}")
```

After saving the file, run in the terminal with `python gene_filter.py`.

### Example 2: Reading a Sample Metadata Table and Filtering Rows 
We’ll now work with a simplified sample_metadata.tsv file, formatted like a real metadata table:

**Goal**: Read a file with tabular data and print only the samples with yield values > 15.0.

Open a file with `vi filter_samples.py`, copy and paste:

```python
# filter_samples.py

# Read the file
with open("sample_metadata.tsv", "r") as infile:
    lines = infile.readlines()

header = lines[0].strip().split("\t")  # Extract column names
data_lines = lines[1:]  # Skip header

print("Samples with treatment = Treated and yield > 15.0:\n")

for line in data_lines:
    parts = line.strip().split("\t")
    sample_id, treatment, yield_val = parts[0], parts[1], float(parts[2])

    if treatment == "Treated" and yield_val > 15.0:
        print(f"{sample_id}: Yield = {yield_val}")
```
After saving the file, run in the terminal with `python filter_samples.py`.

## R and Rstudio

### Open RStudio on the Cluster
On the LiSC cluster, you do not run RStudio directly from the terminal.
Instead, you:

Open a browser and go to:
https://rstudio.lisc.univie.ac.at

Log in with your cluster username and password.
This gives you an RStudio session running on the cluster.

Once logged in, open a new R script and paste:
```r
# ================================================
# INTRO TO R FOR BIOLOGY - DEMO SCRIPT 
# ================================================

# 0. Setup ------------------------------------------------------------
rm(list = ls())

# Install once (if not installed):
# install.packages("ggplot2")

library(ggplot2)

# ------------------------------------------------
# 1. VARIABLES AND DATA TYPES
# ------------------------------------------------
gene_count <- 1500                # numeric
gc_content <- 0.42                # numeric (double)
species <- "Arabidopsis thaliana" # character
is_model <- TRUE                  # logical

class(gene_count)
class(species)

# Vector
chromosomes <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chromosomes

# ------------------------------------------------
# 2. DATA STRUCTURES
# ------------------------------------------------
expr_values <- c(12.5, 15.3, 14.8, 10.1, 9.5)
names(expr_values) <- chromosomes

tissue <- factor(c("Leaf", "Root", "Leaf", "Stem", "Root"))

genes <- data.frame(
  GeneID = paste0("Gene", 1:5),
  Tissue = tissue,
  Expression = expr_values
)
genes

# ------------------------------------------------
# 3. FUNCTIONS
# ------------------------------------------------
gc_content_calc <- function(sequence) {
  g <- sum(strsplit(sequence, "")[[1]] == "G")
  c <- sum(strsplit(sequence, "")[[1]] == "C")
  total <- nchar(sequence)
  return((g + c) / total)
}

gc_content_calc("ATGCGCGTAT")

# ------------------------------------------------
# 4. CONTROL FLOW
# ------------------------------------------------
for (i in 1:nrow(genes)) {
  cat("Gene", genes$GeneID[i], "in", as.character(genes$Tissue[i]),
      "has expression", genes$Expression[i], "\n")
}

if (mean(genes$Expression) > 10) {
  print("Average expression is high")
} else {
  print("Average expression is low")
}

# ------------------------------------------------
# 5. FILE INPUT/OUTPUT
# ------------------------------------------------
write.csv(genes, "genes_demo.csv", row.names = FALSE)
genes_in <- read.csv("genes_demo.csv")
head(genes_in)

# ------------------------------------------------
# 6. PLOTTING WITH ggplot2
# ------------------------------------------------

# Histogram
ggplot(genes, aes(x = Expression)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 5) +
  labs(title = "Histogram of Gene Expression", x = "Expression", y = "Count") +
  theme_minimal()

# Boxplot by tissue
ggplot(genes, aes(x = Tissue, y = Expression, fill = Tissue)) +
  geom_boxplot() +
  labs(title = "Expression by Tissue") +
  theme_minimal()

# Scatterplot
ggplot(genes, aes(x = GeneID, y = Expression, color = Tissue)) +
  geom_point(size = 4) +
  labs(title = "Expression per Gene", x = "Gene", y = "Expression") +
  theme_minimal()

# ------------------------------------------------
# END OF SCRIPT
# ------------------------------------------------
```


## You have completed **Day 2**!

### Useful Scripting Tutorials

- [DataCamp: Bash script tutorial](https://www.datacamp.com/tutorial/how-to-write-bash-script-tutorial) 
- [Geeksforgeeks: Bash scripting](https://www.geeksforgeeks.org/linux-unix/bash-scripting-introduction-to-bash-and-bash-scripting/)
- [Python for Biologists](https://pythonforbiologists.com/tutorial.html)
- [Biopython Tutorial](https://biopython.org/docs/latest/Tutorial/index.html)

