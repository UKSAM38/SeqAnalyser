# SeqAnalyse - Advanced Bioinformatics Tool

A comprehensive bioinformatics application for analyzing sequence files, performing advanced sequence analysis, and generating publication-ready results.

## Features

### Core Analysis Features
- **Multi-format Support**: Load AB1, FASTA, FASTQ, GenBank, EMBL, and text files
- **Advanced Sequence Analysis**: Nucleotide and protein sequence alignment with quality filtering
- **Quality Control**: Quality score visualization and filtering
- **Statistical Analysis**: Comprehensive sequence statistics including GC content, molecular weight, and base composition

### Advanced Features
- **Phylogenetic Analysis**: Create phylogenetic trees using UPGMA method
- **Primer Design**: Automated PCR primer design with quality scoring
- **ORF Analysis**: Open Reading Frame detection and analysis
- **Restriction Analysis**: Restriction enzyme site mapping
- **Multi-frame Translation**: Translation in all 6 reading frames
- **Parallel Processing**: Multi-threaded analysis for improved performance

### Visualization
- **Quality Score Plots**: Visualize sequence quality along the read
- **GC Content Analysis**: GC content distribution plots
- **Amino Acid Composition**: Protein composition analysis
- **Phylogenetic Trees**: Interactive tree visualization
- **Alignment Overview**: Visual alignment summaries

### Export Options
- **PDF Export**: Professional PDF reports with advanced formatting
- **Word Documents**: Formatted Word documents with proper alignment display
- **HTML Export**: Web-compatible HTML reports
- **Multiple Themes**: Dark and light theme support

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

2. **Load Files**:
   - Use "Load Sequence Files" to select your sequence files (supports multiple formats)
   - Use "Load Reference" to select your reference sequence

3. **Configure Analysis**:
   - Set quality thresholds
   - Enable additional analyses (statistics, ORF analysis, restriction sites)
   - Choose analysis type (nucleotide or protein)

4. **Run Analysis**:
   - Click "Analyze Nucleotide" or "Analyze Protein"
   - View results in multiple tabs (Results, Visualization, Phylogenetic, Statistics)

5. **Advanced Features**:
   - Create phylogenetic trees
   - Design PCR primers
   - Export results in multiple formats

## Supported File Formats

### Input Formats
- **AB1**: Chromatogram files from sequencing machines
- **FASTA**: Standard sequence format (.fasta, .fa, .fas)
- **FASTQ**: Sequence format with quality scores (.fastq, .fq)
- **GenBank**: NCBI GenBank format (.gb, .gbk)
- **EMBL**: European Molecular Biology Laboratory format (.embl)
- **Text**: Plain text sequence files (.txt)

### Export Formats
- **PDF**: Professional reports with charts and tables
- **Word**: Microsoft Word documents (.docx)
- **HTML**: Web-compatible reports

## Advanced Analysis Options

### Quality Control
- Configurable quality score thresholds
- Quality-based sequence trimming
- Quality score visualization

### Statistical Analysis
- Base composition analysis
- GC content calculation
- Molecular weight determination
- Dinucleotide frequency analysis

### Phylogenetic Analysis
- Multiple sequence alignment
- Distance matrix calculation
- UPGMA tree construction
- Interactive tree visualization

### Primer Design
- Automated primer pair selection
- Tm and GC content optimization
- Primer dimer checking
- Quality scoring system

## Requirements

- Python 3.8+
- Biopython
- PyQt6
- NumPy
- Matplotlib
- Seaborn
- SciPy
- Scikit-learn
- Python-docx
- FPDF2

## Development

This project is developed using:
- Python for core functionality
- PyQt6 for the GUI
- Biopython for sequence analysis
- Matplotlib/Seaborn for visualization
- Multi-threading for performance optimization

## Performance Features

- **Parallel Processing**: Multi-threaded sequence alignment
- **Memory Optimization**: Efficient handling of large sequence files
- **Quality Filtering**: Smart quality-based sequence processing
- **Batch Processing**: Handle multiple files simultaneously

## User Interface

- **Modern Design**: Dark theme with professional styling
- **Tabbed Interface**: Organized results in multiple tabs
- **Interactive Visualizations**: Matplotlib-based charts and plots
- **Progress Tracking**: Real-time analysis progress
- **Responsive Layout**: Adaptive interface for different screen sizes

## Citation

If you use SeqAnalyse in your research, please cite:
```
SeqAnalyse: Advanced Bioinformatics Tool for Sequence Analysis
Version 2.0, 2025
```
