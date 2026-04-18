# Quinone-Binding Cytochrome Design Analysis

A comprehensive analysis pipeline for protein design, parsing LigandMPNN `.fa` files with confidence scores and Boltz-2 output files. Generates publication-ready figures and statistical summaries for multiple design rounds.

## Features

- **Multi-round analysis**: Process multiple design rounds (`round_1`, `round_2`, etc.) with automatic file discovery
- **LigandMPNN metrics**: Extract sequence confidence, ligand confidence, and sequence recovery from `.fa` files
- **Boltz-2 metrics**: Parse confidence JSON, affinity JSON, and PAE npz files for comprehensive structure quality assessment:
  - pLDDT (per-residue and mean)
  - pTM and ipTM scores
  - Ligand and protein ipTM
  - Complex-level metrics (ipLDDT, PDE, iPDE)
  - Chain-specific pTM and pairwise ipTM
  - Pocket-specific PAE (when binding residues provided)
  - Binding affinity scores
- **Publication-quality figures**: Generates TIFF (600 DPI) and PDF figures with proper formatting for scientific publications
- **Statistical testing**: Automatic t-test comparisons between design rounds with significance annotation
- **Mutation heatmap**: Visualize sequence conservation across rounds relative to round 0
- **Data export**: Comprehensive CSV export of all metrics for downstream analysis

## Installation

### Prerequisites

- Python 3.8 or higher
- Required packages (install via pip):

```bash
pip install numpy scipy matplotlib
```

## Usage

### Command Line Interface

The main analysis script is `main.py`. Run it with the following options:

```bash
python main.py --data-folder /path/to/design/data --output-prefix output_name [--binding-residues RES1 RES2 ...]
```

#### Arguments

- `--data-folder`: Path to the base folder containing design round subfolders (e.g., `round_1`, `round_2`, etc.)
- `--output-prefix`: Prefix for output files (default: `cytochrome_design`)
- `--binding-residues`: Optional list of quinone-binding residue numbers for pocket-specific PAE analysis

### Data Organization

Organize your design data in the following structure:

```
data_folder/
├── round_1/
│   └── ligandmpnn/
│       ├── sample_1.fa
│       ├── sample_2.fa
│       └── ...
├── round_2/
│   └── ligandmpnn/
│       ├── sample_1.fa
│       └── ...
└── round_3/
    └── ligandmpnn/
        └── ...
```

For Boltz-2 data, place the output files in the same round folders:

```
round_1/
├── ligandmpnn/
│   └── ...
├── boltz/
│   ├── confidence.json
│   ├── affinity.json
│   └── pae.npz
└── ...
```

### Output Files

The script generates the following output files:

- **Figures** (TIFF and PDF formats):
  - `{prefix}_ligandmpnn_confidence.tiff/pdf`: LigandMPNN overall confidence
  - `{prefix}_seq_rec.tiff/pdf`: Sequence recovery
  - `{prefix}_boltz_confidence.tiff/pdf`: Boltz-2 confidence score
  - `{prefix}_boltz_plddt.tiff/pdf`: Boltz-2 pLDDT
  - `{prefix}_boltz_ptm_iptm.tiff/pdf`: pTM and ipTM scores
  - `{prefix}_pair_chain_iptm.tiff/pdf`: Pair chain ipTM comparison
  - `{prefix}_mutation_frequency_heatmap.tiff/pdf`: Mutation frequency heatmap
  - Additional figures for affinity and pocket PAE if data is available

- **Data export**:
  - `{prefix}_all_data.csv`: Comprehensive CSV with all metrics

- **Statistics**: Printed to console for figure captions

## Examples

### Basic Analysis

```bash
python main.py --data-folder ./design_rounds --output-prefix my_analysis
```

### With Binding Residues

```bash
python main.py --data-folder ./design_rounds --output-prefix my_analysis --binding-residues 45 67 89 123
```

This enables pocket-specific PAE analysis for the specified residue positions.

## Project Structure

```
statistical_analysis/
├── main.py                 # Main analysis script
├── config.py               # Matplotlib configuration
├── data/
│   ├── __init__.py
│   └── loader.py           # Data loading functions
├── parsers/
│   ├── __init__.py
│   ├── ligandmpnn.py       # LigandMPNN .fa parser
│   └── boltz.py            # Boltz-2 output parsers
├── plots/
│   ├── __init__.py
│   ├── core.py             # Core plotting functions
│   ├── individual.py       # Individual metric plots
│   ├── pair_chain.py       # Pair chain ipTM plots
│   └── heatmap.py          # Mutation heatmap
├── stats/
│   ├── __init__.py
│   └── significance.py     # Statistical significance utilities
└── export/
    ├── __init__.py
    └── csv_export.py       # CSV export functions
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
