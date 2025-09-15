## Week 1: Foundations of Bioinformatics and Genomic Data Handling

## Day 1: Linux, Command Line, and HPC Usage

This session will introduce you to the Linux operating system, which is the backbone of most bioinformatics workflows. We will cover the fundamental concepts of Linux, its file system, and how to interact with it using the command line interface (CLI). A strong understanding of the command line is crucial for efficient data manipulation and running bioinformatics tools.

**What is Linux and why should we use it?**

Linux is a powerful, open-source operating system widely used in scientific computing, including bioinformatics. Its stability, flexibility, and robust command-line tools make it an ideal environment for handling large biological datasets and running computationally intensive analyses.

Linux is the preferred operating system in bioinformatics primarily because it excels at processing and managing enormous volumes of genomic data, ranging from tens of gigabytes to multiple terabytes, containing millions or even billions of sequencing reads. 

Bioinformatics tools such as aligners, SNP callers, RNA‑seq pipelines, and variant callers are designed and optimized specifically for Linux/Unix environments, and many simply won’t build or function correctly on other operating systems. This native compatibility ensures reliable performance when analyzing large-scale datasets.

**Why the Command Line?**

The command line allows you to interact with your computer by typing commands instead of using a graphical user interface (GUI). While it might seem daunting at first, the CLI offers unparalleled control, efficiency, and automation capabilities, which are essential for bioinformatics tasks. Many bioinformatics tools are exclusively command-line based.

**Basic Command Line Operations:**

We will start with the most essential commands that you will use daily:

*   `pwd`: Print Working Directory - Shows your current location in the file system.
*   `ls`: List - Lists the contents of a directory.
    *   `ls -l`: Long format, showing permissions, ownership, size, and modification date.
    *   `ls -a`: All, including hidden files.
*   `cd`: Change Directory - Navigates between directories.
    *   `cd ..`: Move up one directory.
    *   `cd ~`: Go to your home directory.
    *   `cd /path/to/directory`: Go to a specific absolute path.
*   `mkdir`: Make Directory - Creates a new directory.
*   `rmdir`: Remove Directory - Deletes an empty directory.
*   `touch`: Creates an empty file or updates the timestamp of an existing file.
*   `cp`: Copy - Copies files or directories.
    *   `cp file1 file2`: Copies `file1` to `file2`.
    *   `cp -r dir1 dir2`: Recursively copies `dir1` to `dir2`.
*   `mv`: Move/Rename - Moves or renames files or directories.
*   `rm`: Remove - Deletes files or directories.
    *   `rm file`: Deletes a file.
    *   `rm -r dir`: Recursively deletes a directory and its contents.
    *   `rm -f file`: Force delete (no prompt).
    *   `rm -rf dir`: Force recursive delete (use with extreme caution!).
*   `cat`: Concatenate and display files - Displays the content of a file.
*   `less`: Pager - Views file content one screen at a time (useful for large files).
*   `head`: Displays the first few lines of a file.
*   `tail`: Displays the last few lines of a file.
*   `grep`: Global Regular Expression Print - Searches for patterns in text files.
    *   `grep "pattern" filename`: Searches for "pattern" in `filename`.
    *   `grep -i "pattern" filename`: Case-insensitive search.
    *   `grep -r "pattern" directory`: Recursive search in a directory.
*   `|`: Pipe - Redirects the output of one command as the input to another.
*   `>`: Redirects output to a file (overwrites).
*   `>>`: Appends output to a file.

---
## Afternoon Session: Practice Session on Linux

This hands-on session will allow you to apply the basic Linux commands learned in the morning.  
We will work step by step through practical exercises to strengthen your understanding of:

- **File system navigation**
- **File and directory creation**
- **Moving between paths**

The goal is to become comfortable with the Linux environment, since it is the foundation for all subsequent bioinformatics analyses.


**Understanding the Command Prompt**  
When using the terminal, you’ll see a **command prompt** that looks something like:

```bash
[username@server ~]$
```
Where:
- `username`: your login name
- `server`: the name of the machine you're logged into
- `~`: shorthand for your home directory
- `$`: indicates the terminal is ready for your input

**Home Directory: `~` vs `$HOME`**  
To go to your home directory, you can use either:

`cd ~`

or

`cd $HOME`

Both commands work the same. `$HOME` is an environment variable that always stores the path to your home directory. This becomes very useful when writing portable scripts.


### Exercise 1: Navigating the File System

In this first exercise, you will create a working environment for the training and learn how to move between directories.

**Step 1 — Check your current location**  
Always begin by checking where you are in the filesystem with `pwd` (print working directory). This tells you your starting point.  
```bash
pwd
```

**Step 2 — Create a new directory**  
Create a folder called 'bioinformatics_training' inside your home directory (~).
```bash
mkdir ~/bioinformatics_training
```

**Step 3 — Enter the new directory**  
Move into the new folder using cd (change directory).
```bash
cd ~/bioinformatics_training
```

**Step 4 — Create subdirectories**  
Inside this folder, create two subfolders:  
`data`    → where we will keep files and datasets  
`scripts` → where we will write small shell scripts later

```bash
mkdir data scripts
```

**Step 5 — Navigate into a subdirectory**  
Enter the 'data' directory.
```bash
cd data
```

**Step 6 — Return to the parent directory**  
Use cd .. to go one level up and return to 'bioinformatics_training'.
```bash
cd ..
```

**Step 7 — Confirm your location**  
Check where you are in the filesystem with pwd (print working directory).
```bash
pwd
```

### Exercise 2: Creating and Manipulating Files

In this exercise, you will learn how to **create**, **edit**, **view**, **copy**, **move**, and **rename** files.  
These operations are the foundation of working with any dataset in Linux, whether you are handling sequencing reads, annotation files, or analysis outputs.

We will work inside the folder `~/bioinformatics_training/data` created in Exercise 1.


**Step 1 — Go to the data directory**  
We first move into the `data` folder where we want to work.  
The command `cd` means "change directory".
```bash
cd ~/bioinformatics_training/data
```

Always check where you are with `pwd` (print working directory). *This helps avoid mistakes like creating files in the wrong folder.
```bash
pwd
```
Expected: `/home/<username>/bioinformatics_training/data`

**Step 2 — Create an empty file**  
The command `touch` creates an empty file if it does not exist.
Here we create a file called sample_data.txt.
```bash
touch sample_data.txt
```

Use `ls -l` to check that the file was created.
The file will have size 0 since it is empty.
```bash
ls -l sample_data.txt
```

**Step 3 — Add text to the file**  
The `echo` command prints text. Using `>` sends the text into a file.
If the file already exists, `>` will overwrite its contents.
```bash
echo "This is the first line of sample data." > sample_data.txt
```

To add more lines without overwriting, use `>>` (append).
```bash
echo "This is the second line." >> sample_data.txt
```

**Step 4 — View the file contents**  
`cat` shows the whole file on the screen (good for small files).
```bash
cat sample_data.txt
```

For larger files, `less` is safer: it lets you scroll up and down.
Press `q` to quit less.
```bash
less sample_data.txt
```

**Step 5 — Copy the file**  
The command `cp` makes a copy of a file. Here we copy sample_data.txt
into a new file called data_backup.txt.
```bash
cp sample_data.txt data_backup.txt
```

Verify that both files exist and compare their sizes.
```bash
ls -l sample_data.txt data_backup.txt
```

**Step 6 — Move the copy to the scripts directory**  
The command `mv` can move a file from one location to another.  
`..` means "parent directory". So `../scripts/` is the scripts folder
one level above the current folder.
```bash
mv data_backup.txt ../scripts/
```

Check that the file has indeed moved to scripts.
```bash
ls -l ../scripts/
```

**Step 7 — Rename inside scripts**  
The same `mv` command is also used to rename files.
Here we rename data_backup.txt to first_script_notes.txt.
```bash
cd ../scripts/
mv data_backup.txt first_script_notes.txt
```

Verify the rename worked.
```bash
ls -l
```

**Step 8 — Verify the final state**  
Finally, check that:
- The original sample_data.txt is still inside data.
- The renamed file first_script_notes.txt is inside scripts.
```bash
cd ../data
ls -l
```

### Exercise 3: Searching and Filtering Text

In this exercise, you will learn how to **search and filter information** inside text files.  
These skills are essential in bioinformatics because almost all genomic data (FASTA, FASTQ, GFF, VCF, etc.) are stored as plain text.

We will practice with the command `grep` (search patterns in text), and combine `head` and `tail` to extract specific parts of files.

**Step 1 — Go to the data directory**  
```bash
cd ~/bioinformatics_training/data
pwd
```
Expected: `/home/<username>/bioinformatics_training/data`

**Step 2 — Create a test file with genomic sequences**  
We will create a small file called genomic_sequences.txt that looks like a FASTA file.  
Lines starting with `>` are sequence identifiers, followed by DNA sequences.
```bash
echo -e ">gene1\nATGCATGCATGC\n>gene2\nCGTACGTACGTAG\nATGCATGC" > genomic_sequences.txt
```

Check that the file exists and view it:
```bash
ls -l genomic_sequences.txt
cat genomic_sequences.txt
```

**Step 3 — Search for sequence identifiers**  
`grep ">" file` finds all lines that contain `>`.  
In FASTA files, identifiers always begin with `>`.
```bash
grep ">" genomic_sequences.txt
```

**Step 4 — Search for a DNA motif**  
Here we search for lines that contain the motif `ATGC`. 
This simulates searching for specific sequence patterns.
```bash
grep "ATGC" genomic_sequences.txt
```
Expected output will show lines where ATGC appears.

**Step 5 — Extract a specific line using head and tail**  
`head -n N` prints the first N lines of a file.  
`tail -n 1` takes the last line of whatever head printed.
Combined, we can extract line 1, line 2, etc:
```bash
head -n 1 genomic_sequences.txt | tail -n 1
```

Example: get the 3rd line:
```bash
head -n 3 genomic_sequences.txt | tail -n 1
```
Expected: the third line ("ATGCATGCATGC")


### Exercise 4: Working with Compressed Files and Downloading Data

#### Step 1 — Navigate to your `data` folder:
```bash
cd ~/bioinformatics_training/data
```

#### Step 2 — Download a sample file:
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/README.txt
# OR using curl
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/README.txt
```

You can inspect the file with
```bash
cat README.txt
```

#### Step 3 — Download a compressed FASTA file:
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
```

#### Step 4 — Inspect the `.gz` file:
```bash
ls -lh *.gz
zcat GCF_000001735.4_TAIR10.1_genomic.fna.gz | head
# OR
gunzip -c GCF_000001735.4_TAIR10.1_genomic.fna.gz | head
```

#### Step 5 — Decompress the file:
```bash
gunzip -c GCF_000001735.4_TAIR10.1_genomic.fna.gz > genome.fna
```

#### Step 6 — Verify the decompressed file:
```bash
ls -lh genome.fna
head genome.fna
```

---
## Introduction to HPC Usage 

While direct hands-on HPC (High-Performance Computing) usage will be covered in more detail later (e.g., Day 5 with Apocrita), it's important to understand the concept from Day 1. HPC clusters are powerful computing resources used for large-scale bioinformatics analyses that require significant computational power and memory. They typically involve:

  - **Login Nodes:** Where you log in and prepare your jobs.  
  - **Compute Nodes:** Where your actual analyses run.  
  - **Job Schedulers (e.g., Slurm, PBS):** Systems that manage and allocate resources for your jobs.  
  - **Shared File Systems:** Where data is stored and accessed by all nodes.  

For this course, we will simulate some HPC interactions using the command line, and later, if applicable, we will connect to a real HPC environment like Apocrita. The commands you learn today are directly transferable to an HPC environment.

   **Key HPC Concepts:**  
  - **Batch Jobs:** Submitting scripts to run on compute nodes without direct interaction.  
  - **Resource Allocation:** Requesting specific CPU cores, memory, and time for your jobs.  
  - **Modules:** Software management system to load different versions of bioinformatics tools.  

We will delve deeper into these concepts and practical applications in later sessions, especially when we discuss data quality control and large-scale analyses. For now, focus on mastering the basic Linux commands, as they are the foundation for interacting with any computing environment, including HPCs.


### Exercise 4: Introduction to HPC Usage (LISC system, UniVie)

High-Performance Computing (HPC) clusters like **LISC** at the University of Vienna are used for large-scale bioinformatics analyses.  
Instead of running heavy computations on your laptop, you send jobs to the cluster, which has many powerful compute nodes.

In this exercise, you will learn the **basic commands** to interact with LISC using the **Slurm workload manager**.


**Step 1 — Connect to the HPC system**  
Normally you log in via SSH from your terminal.
```bash
ssh <username>@login02.lisc.univie.ac.at
```

After login, check where you are:
```bash
pwd
```
Expected: /lisc/scratch/<username> or similar

**Step 2 — Check basic information**  
Show available resources (queues/partitions):
```bash
sinfo
```

Show your current jobs (none yet):
```bash
squeue -u <username>
```

**Step 3 — Load software with modules**  
HPC systems don’t have all tools loaded by default.
Instead, they use "modules" to load specific software.
```bash
module avail     # list available modules
module load samtools/1.15
samtools --version
module unload samtools/1.15
```

## You have completed **Day 1**!

### Useful Linux Command Tutorials

- [Ubuntu: Command Line for Beginners](https://ubuntu.com/tutorials/command-line-for-beginners#1-overview)  
- [DigitalOcean: Linux Commands](https://www.digitalocean.com/community/tutorials/linux-commands)  
- [Hostinger: Linux Commands](https://www.hostinger.com/tutorials/linux-commands)
- [Shell Genomics](https://datacarpentry.github.io/shell-genomics/)  
