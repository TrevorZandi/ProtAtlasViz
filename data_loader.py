"""
Data loading and processing module for ProtAtlasViz
Handles loading RNA tissue consensus data and histology dictionary
"""
import pandas as pd
from typing import Dict, List, Tuple


class DataLoader:
    """Loads and processes gene expression and tissue grouping data"""
    
    def __init__(self, data_path: str = 'data/rna_tissue_consensus.tsv', 
                 histology_path: str = 'data/Histology_Dictionary.txt'):
        self.data_path = data_path
        self.histology_path = histology_path
        self.expression_data = None
        self.tissue_groups = None
        self.genes = None
        
    def load_expression_data(self) -> pd.DataFrame:
        """Load RNA tissue consensus data"""
        self.expression_data = pd.read_csv(self.data_path, sep='\t')
        # Get unique genes for autocomplete
        self.genes = sorted(self.expression_data['Gene name'].unique())
        return self.expression_data
    
    def load_histology_dictionary(self) -> Dict[str, List[str]]:
        """Parse histology dictionary to create organ -> tissues mapping"""
        tissue_groups = {}
        current_group = None
        
        with open(self.histology_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                # Lines starting with # are group headers
                if line.startswith('#'):
                    current_group = line[1:].strip()
                    tissue_groups[current_group] = []
                elif current_group:
                    tissue_groups[current_group].append(line)
        
        self.tissue_groups = tissue_groups
        return tissue_groups
    
    def get_tissue_to_group_mapping(self) -> Dict[str, str]:
        """Create reverse mapping: tissue -> organ group (case-insensitive)"""
        if self.tissue_groups is None:
            self.load_histology_dictionary()
        
        if self.expression_data is None:
            self.load_expression_data()
        
        # Get all tissues from actual data
        all_tissues_in_data = self.expression_data['Tissue'].unique()
        
        tissue_to_group = {}
        for group, tissues in self.tissue_groups.items():
            for tissue in tissues:
                # Find matching tissue in data (case-insensitive)
                matching_tissues = [t for t in all_tissues_in_data 
                                   if t.lower() == tissue.lower()]
                if matching_tissues:
                    # Map the actual tissue name from data to the group
                    tissue_to_group[matching_tissues[0]] = group
        
        return tissue_to_group
    
    def get_expression_for_genes(self, gene_names: List[str]) -> pd.DataFrame:
        """Get expression data for selected genes"""
        if self.expression_data is None:
            self.load_expression_data()
        
        filtered_data = self.expression_data[
            self.expression_data['Gene name'].isin(gene_names)
        ].copy()
        
        return filtered_data
    
    def get_ordered_tissues_by_group(self) -> Tuple[List[str], Dict[str, str]]:
        """
        Get tissues ordered by organ group
        Returns: (ordered_tissues, tissue_to_group_dict)
        """
        if self.tissue_groups is None:
            self.load_histology_dictionary()
        
        if self.expression_data is None:
            self.load_expression_data()
        
        tissue_to_group = self.get_tissue_to_group_mapping()
        
        # Get all tissues from data that are in the dictionary
        all_tissues_in_data = self.expression_data['Tissue'].unique()
        
        # Order tissues by group
        ordered_tissues = []
        for group, tissues in self.tissue_groups.items():
            for tissue in tissues:
                # Normalize tissue names (case-insensitive matching)
                matching_tissues = [t for t in all_tissues_in_data 
                                   if t.lower() == tissue.lower()]
                if matching_tissues:
                    ordered_tissues.append(matching_tissues[0])
        
        # Add any tissues not in dictionary at the end
        for tissue in all_tissues_in_data:
            if tissue not in ordered_tissues:
                ordered_tissues.append(tissue)
                tissue_to_group[tissue] = "Other"
        
        return ordered_tissues, tissue_to_group
    
    def aggregate_by_organ_group(self, data: pd.DataFrame) -> pd.DataFrame:
        """Aggregate expression data by organ groups"""
        tissue_to_group = self.get_tissue_to_group_mapping()
        
        # Add organ group column
        data = data.copy()
        data['Organ Group'] = data['Tissue'].map(tissue_to_group)
        data['Organ Group'] = data['Organ Group'].fillna('Other')
        
        # Aggregate by gene and organ group (mean nTPM)
        aggregated = data.groupby(['Gene', 'Gene name', 'Organ Group'])['nTPM'].mean().reset_index()
        aggregated.rename(columns={'Organ Group': 'Tissue'}, inplace=True)
        
        return aggregated
