**Overview**

This repository provides a Python-based automation pipeline for high-throughput analysis of Sanger sequencing chromatogram files (.ab1) using the INDIGO webserver from Gear Genomics.

The pipeline leverages Selenium WebDriver to automate:

1) Uploading .ab1 sample files

2) Uploading wild-type reference files

3) Running INDIGO analysis

4) Extracting HTML result pages

5)Highlighting Protospacer Adjacent Motif (PAM) regions

6) Batch processing multiple gRNA targets

This tool transforms the INDIGO web application into a scalable solution for CRISPR/Cas9 mutation validation in plant genome editing workflows.

**Background**

CRISPR/Cas9 genome editing generates insertions/deletions (indels) that must be validated through sequencing. While Sanger sequencing is cost-effective for screening large plant populations, manual analysis of .ab1 files through the INDIGO web interface is:

Time-consuming

Repetitive

Error-prone

Difficult to scale

This automation pipeline eliminates manual intervention and enables high-throughput mutation analysis.

**Key Features**

Batch processing of multiple .ab1 files

Automated upload to INDIGO webserver

Wild-type reference comparison

PAM sequence highlighting in HTML reports

Support for multiple gRNA regions

Error handling for missing files and web failures

Reproducible and scalable workflow

**Workflow**

Input Validation

Check for existence of input folder containing .ab1 files

File Upload

Upload sample and wild-type reference files to INDIGO

Automated Analysis

Extract page source

Identify and highlight PAM regions

Result Retrieval

Save INDIGO-generated HTML reports

Cleanup

Close browser session

Repeat

Process all files in the input directory

**System Requirements**

Tested on:

Windows 10

16 GB RAM

Intel® Core™ i5 processor

Google Chrome (v126.0)

Compatible ChromeDriver

Python 3.8

Dependencies

**Install required Python packages:**

pip install selenium

**You must also install:**

Google Chrome

Matching ChromeDriver

Python 3.8+

Directory Structure
project_root/
│
├── input_files/          # Edited sample .ab1 files
├── wild_type/            # Wild-type reference .ab1 file
├── output/               # Generated HTML reports
├── script.py             # Automation script
└── README.md
Usage

Place edited .ab1 files inside the input_files/ directory

Place wild-type reference .ab1 file inside the wild_type/ directory

Update PAM region coordinates inside the script (if needed)

**Run:**

python script.py

Processed HTML reports will be saved in the output folder.

Performance
Accuracy

100% concordance with manual INDIGO analysis

Accurate PAM region highlighting

Reliable indel detection

Efficiency

Processes large numbers of .ab1 files automatically

Eliminates manual upload and inspection

Handles multiple gRNAs in a single run

Reliability

Built-in exception handling

Logs errors without halting entire workflow

Reproducible results across runs

Example Output

INDIGO alignment reports (Alt1 & Alt2)

PAM regions highlighted

gRNA-specific indels clearly visible

HTML output files ready for documentation or publication

**Applications**

CRISPR/Cas9 mutation validation

High-throughput plant genome editing workflows

Large-scale screening of edited crop lines

Functional genomics studies

Future Improvements

Parallel processing support

Additional sequencing format support

Integration with other mutation analysis tools

Automated report summarization

**Citation**

If you use this pipeline in your research, please cite:

Suresh, V., Girish, C., & Tavva, V.S.S.
Python-based automation of INDIGO webserver using Selenium: A high throughput analysis of Sanger sequence data to detect allelic variations created by CRISPR/Cas9-mediated genome editing of crop plants.
bioRxiv (2025).

