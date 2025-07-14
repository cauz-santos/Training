## Week 1: Foundations of Bioinformatics and Genomic Data Handling

### Day 1: Linux, Command Line, and HPC Usage

This session will introduce you to the Linux operating system, which is the backbone of most bioinformatics workflows. We will cover the fundamental concepts of Linux, its file system, and how to interact with it using the command line interface (CLI). A strong understanding of the command line is crucial for efficient data manipulation and running bioinformatics tools.

**What is Linux and why should we use it?**

Linux is a powerful, open-source operating system widely used in scientific computing, including bioinformatics. Its stability, flexibility, and robust command-line tools make it an ideal environment for handling large biological datasets and running computationally intensive analyses.
Linux is the preferred operating system in bioinformatics primarily because it excels at processing and managing enormous volumes of genomic data—ranging from tens of gigabytes to multiple terabytes, containing millions or even billions of sequencing reads. Bioinformatics tools such as aligners, SNP callers, RNA‑seq pipelines, and variant callers are designed and optimized specifically for Linux/Unix environments, and many simply won’t build or function correctly on other operating systems. This native compatibility ensures reliable performance when analyzing large-scale datasets.

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


### Day 1: Linux, Command Line, and HPC Usage (Continued)

**Afternoon Session: Practice Session on Linux**

This hands-on session will allow you to apply the basic Linux commands learned in the morning. We will work through a series of exercises designed to solidify your understanding of file system navigation, file manipulation, and basic text processing. The goal is to build your confidence and proficiency in using the command line, which is fundamental for all subsequent bioinformatics tasks.

**Practice Exercises:**

1.  **Navigating the File System:**
    *   Create a new directory named `bioinformatics_training` in your home directory.
    *   Navigate into this new directory.
    *   Inside `bioinformatics_training`, create two subdirectories: `data` and `scripts`.
    *   Navigate into the `data` directory.
    *   Return to the `bioinformatics_training` directory.

    ```bash
    mkdir ~/bioinformatics_training
    cd ~/bioinformatics_training
    mkdir data scripts
    cd data
    cd ..
    pwd
    ```

2.  **Creating and Manipulating Files:**
    *   Inside the `data` directory, create an empty file named `sample_data.txt`.
    *   Add some text to `sample_data.txt` using the `echo` command and redirection (e.g., `echo "This is a test line." > sample_data.txt`).
    *   View the content of `sample_data.txt` using `cat` and `less`.
    *   Copy `sample_data.txt` to `data_backup.txt` within the same directory.
    *   Move `data_backup.txt` from the `data` directory to the `scripts` directory.
    *   Rename `data_backup.txt` in the `scripts` directory to `first_script_notes.txt`.

    ```bash
    cd data
    touch sample_data.txt
    echo "This is the first line of sample data." > sample_data.txt
    echo "This is the second line." >> sample_data.txt
    cat sample_data.txt
    less sample_data.txt
    cp sample_data.txt data_backup.txt
    mv data_backup.txt ../scripts/
    cd ../scripts/
    mv data_backup.txt first_script_notes.txt
    ls -l
    ```

3.  **Searching and Filtering Text:**
    *   Create a file named `genomic_sequences.txt` in the `data` directory with some example sequences (e.g., lines starting with `>gene1`, `>gene2`, `ATGCATGC`, `CGTACGTA`).
    *   Use `grep` to find all lines that start with `>`.
    *   Use `grep` to find all lines containing `ATGC`.
    *   Combine `head` and `tail` with pipes (`|`) to view specific parts of a file.

    ```bash
    cd ../data
    echo ">gene1\nATGCATGCATGC\n>gene2\nCGTACGTACGTAG\nATGCATGC" > genomic_sequences.txt
    grep ">" genomic_sequences.txt
    grep "ATGC" genomic_sequences.txt
    head -n 1 genomic_sequences.txt | tail -n 1
    ```

4.  **Introduction to HPC Usage (Conceptual):**

    While direct hands-on HPC (High-Performance Computing) usage will be covered in more detail later (e.g., Day 5 with Apocrita), it's important to understand the concept from Day 1. HPC clusters are powerful computing resources used for large-scale bioinformatics analyses that require significant computational power and memory. They typically involve:

    *   **Login Nodes:** Where you log in and prepare your jobs.
    *   **Compute Nodes:** Where your actual analyses run.
    *   **Job Schedulers (e.g., Slurm, PBS):** Systems that manage and allocate resources for your jobs.
    *   **Shared File Systems:** Where data is stored and accessed by all nodes.

    For this course, we will simulate some HPC interactions using the command line, and later, if applicable, we will connect to a real HPC environment like Apocrita. The commands you learn today are directly transferable to an HPC environment.

    **Key HPC Concepts:**
    *   **Batch Jobs:** Submitting scripts to run on compute nodes without direct interaction.
    *   **Resource Allocation:** Requesting specific CPU cores, memory, and time for your jobs.
    *   **Modules:** Software management system to load different versions of bioinformatics tools.

    We will delve deeper into these concepts and practical applications in later sessions, especially when we discuss data quality control and large-scale analyses. For now, focus on mastering the basic Linux commands, as they are the foundation for interacting with any computing environment, including HPCs.

