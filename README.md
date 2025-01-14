# Liquid-DNA-Drive-System
# All MATLAB file

This repository contains a MATLAB script for processing DNA sequences from FASTQ files. It performs sequence assembly, quality score analysis, and includes functionalities for handling oligonucleotide data.

---

## 1. System Requirements

### 1.1 Software Dependencies and Operating Systems
- **MATLAB**: R2021b or R2023a
- **Bioinformatics Toolbox** for MATLAB (required for FASTQ file processing)

### 1.2 Versions the Software Has Been Tested On
- MATLAB R2021b on Windows 10 Pro
- MATLAB R2021b on macOS Sonoma 14.7

### 1.3 Required Non-Standard Hardware
- Standard desktop or laptop with minimum 8 GB RAM.
- No additional hardware is required. 

---

## 2. Installation Guide

### 2.1 Instructions
1. Clone this repository to your local machine:
git clone https://github.com/bioelectronicsbiooptics/Liquid-DNA-Drive-System/blob/main
2. Ensure MATLAB and Bioinformatics Toolbox are installed.
3. Open the script file fastq_analysis_2.m in MATLAB.

### 2.2 Typical Install Time
- Approximately 2-3 minutes on a standard desktop computer.

---

## 3. Demo

###3.1 Instructions to Run on Data
1. Prepare two FASTQ files (e.g., R1.fastq and R2.fastq) in the same directory as the script.
2. Modify the script to include your file names:
[~, S1, Q1] = fastqread('R1.fastq');
[~, S2, Q2] = fastqread('R2.fastq'); 
3. Run the script in MATLAB.

### 3.2 Expected Output
- Assembled sequences stored in variables S1S and S1Q.
- Visualization of the histogram graph for sequences with length trimming.
- Visualization the nucleotide proportion at each position with primer matching, quality trimming, length trimming.

### 3.3 Expected fun time for demo on a “normal” desktop computer
- Approximately 10-20 minutes for FASTQ files on a standard desktop computer.

---

## 4. Instruction for Use

### 4.1 How to Run the Software on Your Data
1. Place your FASTQ files in the working directory.
2. Update the script with your FASTQ file manes.
3. Run the script in MATLAB, following any prompts (e.g., inputting sequence lengths).
