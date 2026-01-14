# 🧬 ProtAtlasViz

An interactive web application for visualizing tissue-selective gene expression patterns from RNA-seq data. Built with Dash and Plotly, this tool allows researchers to explore expression profiles of multiple genes across different tissues and organ systems.

## Features

- **Multi-Gene Selection**: Select up to 10 genes using an autocomplete search interface
- **Multiple Visualizations**: 
  - Heatmap (default) - Shows expression density across all tissues
  - Grouped Bar Chart - Compare expression levels across genes
  - Box Plot - View expression distributions
- **Hierarchical Tissue Organization**: 
  - View individual tissues organized by organ systems
  - Aggregate data at the organ group level
- **Interactive**: Hover over data points for detailed information
- **Responsive Design**: Bootstrap-based UI that works on various screen sizes

## Data Source

This application uses RNA tissue consensus data from The Human Protein Atlas (www.proteinatlas.org), containing normalized transcript expression values (nTPM) across multiple human tissues.

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Setup

1. Clone the repository:
```bash
git clone https://github.com/TrevorZandi/ProtAtlasViz.git
cd ProtAtlasViz
```

2. Create and activate a virtual environment (recommended):
```bash
python -m venv .venv
source .venv/bin/activate  # On Linux/Mac
# or
.venv\Scripts\activate  # On Windows
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

1. Start the application:
```bash
python app.py
```

2. Open your web browser and navigate to:
```
http://127.0.0.1:8050/
```

3. Use the interface:
   - **Gene Selection**: Type to search and select genes (up to 10)
   - **Grouping Level**: Choose between individual tissues or organ groups
   - **Visualization Type**: Select your preferred chart type
   - The visualization updates automatically based on your selections

## Project Structure

```
ProtAtlasViz/
├── app.py                          # Main Dash application
├── data_loader.py                  # Data loading and processing module
├── requirements.txt                # Python dependencies
├── README.md                       # This file
├── data/
│   ├── rna_tissue_consensus.tsv   # RNA expression data
│   └── Histology_Dictionary.txt   # Tissue-to-organ grouping
└── .gitignore                      # Git ignore file
```

## Dependencies

- dash - Web application framework
- plotly - Interactive visualization library
- pandas - Data manipulation and analysis
- numpy - Numerical computing
- dash-bootstrap-components - Bootstrap components for Dash

## Development

To contribute or modify the application:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

## License

Apache License 2.0 - See LICENSE file for details

## Acknowledgments

- Data sourced from The Human Protein Atlas (www.proteinatlas.org)
- Built with Dash by Plotly

## Contact

For questions or issues, please open an issue on GitHub.
