## Week 1: Foundations of Bioinformatics and Genomic Data Handling

### Day 2: Bioinformatics Scripting Essentials (Bash, Python)

**Session: Writing Bioinformatics Scripts (Parsing FASTA Files, Automation)**

On Day 2, we transition from using individual commands to combining them into scripts. Scripting is the key to automating repetitive tasks, creating reproducible workflows, and handling large datasets efficiently. We will cover the fundamentals of scripting in two powerful languages for bioinformatics: Bash and Python.

**Introduction to Bash Scripting**

Bash (Bourne-Again SHell) is the default command-line interpreter on most Linux systems. A Bash script is simply a text file containing a series of commands that are executed sequentially. It's an excellent tool for automating command-line workflows.

**Creating a Simple Bash Script:**

1.  **Create a script file:** Use a text editor like `nano` or `vim` to create a file, for example, `my_first_script.sh`.
2.  **Add the shebang:** The first line of a Bash script should be `#!/bin/bash`. This tells the system to use the Bash interpreter to run the script.
3.  **Add commands:** Add the commands you want to execute, one per line.
4.  **Make the script executable:** Use the `chmod +x my_first_script.sh` command to give the script execute permissions.
5.  **Run the script:** Execute the script by typing `./my_first_script.sh`.

**Example Bash Script for Automation:**

Let's create a script that automates the process of creating a new project directory with subdirectories for data, scripts, and results.

```bash
#!/bin/bash

# This script creates a new project directory structure

PROJECT_NAME=$1 # Get the project name from the command line

mkdir $PROJECT_NAME
cd $PROJECT_NAME
mkdir data scripts results

echo "Project directory 
```

To run this script, you would save it as `create_project.sh`, make it executable (`chmod +x create_project.sh`), and then run it with a project name as an argument: `./create_project.sh my_new_project`.

For example:
```bash
./create_project.sh my_new_project Test
```


**Introduction to Python for Bioinformatics**

Python is a versatile, high-level programming language that has become a cornerstone of bioinformatics due to its readability, extensive libraries (like Biopython), and powerful data analysis capabilities. We will cover the absolute basics to get you started.

**Python Fundamentals:**

*   **Data Types:** We will introduce fundamental data types such as strings (for sequences), integers, floats, lists (for storing collections of items), and dictionaries (for key-value pairs).
*   **Loops:** Loops are essential for iterating over data. We will cover `for` loops, which are used to iterate over a sequence (like a list or a string).
*   **File Handling:** A crucial skill in bioinformatics is reading from and writing to files. We will learn how to open files, read their content line by line, and write results to new files.

**Parsing FASTA Files with Python:**

FASTA is a common text-based format for representing nucleotide or peptide sequences. A FASTA file can contain one or more sequences, each with a header line starting with `>` followed by the sequence data on subsequent lines.

Here is a simple Python script to parse a FASTA file and print the header and sequence for each entry:

```python
# fasta_parser.py

def parse_fasta(filename):
    """Parses a FASTA file and yields header, sequence pairs."""
    with open(filename, 'r') as f:
        header = ''
        sequence = ''
        for line in f:
            line = line.strip() # Remove leading/trailing whitespace
            if line.startswith('>'):
                if header: # If we have a previous sequence, yield it
                    yield header, sequence
                header = line[1:] # Remove the '>'
                sequence = ''
            else:
                sequence += line
        if header: # Yield the last sequence in the file
            yield header, sequence

# Example usage:
if __name__ == '__main__':
    fasta_file = 'example.fasta' # Assume this file exists
    # Create a dummy fasta file for the example
    with open(fasta_file, 'w') as f:
        f.write('>gene1 description of gene1\n')
        f.write('ATGCATGCATGC\n')
        f.write('GATTACAGATTACA\n')
        f.write('>gene2 description of gene2\n')
        f.write('CGTACGTACGTACG\n')

    for header, sequence in parse_fasta(fasta_file):
        print(f'Header: {header}')
        print(f'Sequence: {sequence}')
        print(f'Sequence Length: {len(sequence)}\n')

```

This script demonstrates key Python concepts: functions, file handling, loops, and string manipulation. During the hands-on session, you will write and run this script, and we will explore how to modify it to perform other tasks, such as calculating GC content or searching for specific motifs in the sequences.

By the end of Day 2, you will have a foundational understanding of how to automate bioinformatics workflows using both Bash and Python, with a practical focus on parsing one of the most common data formats in genomics.


