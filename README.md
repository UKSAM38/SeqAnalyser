# SeqAnalyse - Bioinformatics Tool

A Mac application for analyzing AB1 sequence files, performing BLAST analysis, and protein sequence analysis.

## Features

- Load and read AB1 format files
- Perform nucleotide BLAST analysis
- Translate DNA sequences to protein sequences
- Perform protein sequence alignment
- Visualize sequence differences

## Installation

1. Ensure you have Python 3.8 or higher installed
2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

1. Run the application:
```bash
python src/main.py
```

2. Use the "Load AB1 File(s)" button to select your sequence files
3. Choose the analysis type you want to perform:
   - BLAST Analysis
   - Protein Translation
   - Sequence Alignment

## Requirements

- Python 3.8+
- Biopython
- PyQt6
- numpy
- matplotlib

## Development

This project is developed using:
- Python for core functionality
- PyQt6 for the GUI
- Biopython for sequence analysis
