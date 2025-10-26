# Multiple Protein Sequence Alignment App

A Python-based desktop application for multiple protein sequence alignment using MUSCLE and ClustalW algorithms.

## Overview

This application provides an intuitive tool for analyzing protein sequences through multiple sequence alignment (MSA). It integrates MUSCLE and ClustalW algorithms, enabling users to compare alignment quality, identify conserved regions, and explore evolutionary relationships.

## Features

- **Dual Algorithm Support**: MUSCLE and ClustalW
- **User-Friendly Interface**: Streamlit-based web interface
- **Visualizations**: Color-coded alignments, similarity matrices, conservation heatmaps, phylogenetic trees
- **Example Data**: Pre-loaded sequences for quick testing
- **Downloadable Results**: Export alignments, guide trees, and matrices
- **Cross-Platform**: Windows and Linux support

## Quick Start

1. Extract ZIP file to desired location
2. Navigate to `build/exe.win-amd64-3.11/`
3. Run `app.exe`
4. Press Enter if prompted for email
5. Application opens in your browser

## Usage

1. **Upload** your FASTA file or click "Use Example"
2. **Select** algorithm (MUSCLE or ClustalW)
3. **Run Alignment**
4. **View Results** in five tabs:
   - Tool Output: Aligned sequences with color coding
   - Conserved Regions: Similarity matrix and heatmap
   - Dendrogram: Phylogenetic tree
   - Evaluation: Algorithm metrics
   - Download: Export results

## Input Format

FASTA format protein sequences:
```
>SEQUENCE_NAME Description
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
>SEQUENCE_NAME2 Description
MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHG
```

## Conservation Symbols

- `*` = Complete conservation
- `:` = Strong conservation
- `.` = Weak conservation
- ` ` = No conservation

## Technical Stack

- Python, Streamlit, Biopython
- MUSCLE v3.8.31, ClustalW2, FastTree

## Output Files

- Aligned Sequences (.fasta, .aln)
- Guide Trees (.dnd)
- Similarity Matrix (.txt)

## Use Cases

- Evolutionary analysis
- Conserved region identification
- Protein function prediction
- Structural analysis
- Drug discovery applications

## Author

**Shaikh Basim Mushtaque Ahmed**  
Biomedical Engineering - Information Systems in Medicine  
Faculty of Biomedical Engineering, Zabrze, 2024  
Supervisor: Anna Tamulewicz, PhD Eng.

## References

- MUSCLE: https://drive5.com/muscle/
- ClustalW: http://www.clustal.org/clustal2/
- FastTree: http://www.microbesonline.org/fasttree/
- Biopython: https://biopython.org/

---

**Note**: This application is designed for educational and research purposes.Retry
