#!/usr/bin/env python3
"""
Analysis script for quinone-binding membrane cytochrome design
Parses LigandMPNN .fa files with confidence scores and Boltz-2 output files
Includes Boltz confidence metrics
"""

import argparse
from scipy import stats

from config import rcParams  # This applies the matplotlib settings
from data.loader import load_quinone_design_data
from plots.individual import plot_all_metrics
from plots.heatmap import plot_mutation_heatmap
from export.csv_export import export_all_data, print_statistics
from stats.significance import format_significance_label


def main():
    parser = argparse.ArgumentParser(description='Analyze quinone-binding cytochrome design data with all Boltz metrics')
    parser.add_argument('--data-folder', type=str, required=True,
                       help='Base folder containing design round subfolders')
    parser.add_argument('--output-prefix', type=str, default='cytochrome_design',
                       help='Prefix for output files')
    parser.add_argument('--binding-residues', type=int, nargs='+',
                       help='Quinone-binding residue numbers (space-separated)')
    
    args = parser.parse_args()
    
    # Load all design data
    design_data = load_quinone_design_data(
        args.data_folder,
        binding_residues=args.binding_residues
    )
    
    if not design_data['round_names']:
        print("No design data found!")
        return
    
    # Calculate significance annotations for key metrics
    significance_annotations = {}
    round_names = design_data['round_names']
    
    # Only proceed if we have at least 2 rounds with data
    if len(round_names) >= 2:
        # Dictionary to store which metrics have data
        metrics_to_test = {
            'ligandmpnn_confidence': design_data['ligandmpnn_confidence'],
            'ligandmpnn_seq_rec': design_data['ligandmpnn_seq_rec'],
            'boltz_confidence': design_data['boltz_confidence'],
            'boltz_plddt': design_data['boltz_plddt'],
            'boltz_ptm': design_data['boltz_ptm'],
            'boltz_iptm': design_data['boltz_iptm'],
            'boltz_ligand_iptm': design_data['boltz_ligand_iptm'],
            'boltz_protein_iptm': design_data['boltz_protein_iptm'],
            'boltz_complex_iplddt': design_data['boltz_complex_iplddt'],
            'boltz_complex_pde': design_data['boltz_complex_pde'],
            'boltz_complex_ipde': design_data['boltz_complex_ipde'],
            'boltz_chains_ptm_mean': design_data['boltz_chains_ptm_mean'],
            'boltz_pair_iptm_mean': design_data['boltz_pair_iptm_mean'],
            'boltz_pocket_pae': design_data['boltz_pocket_pae'],
            'boltz_mean_pae_npz': design_data['boltz_mean_pae_npz'],
            'boltz_affinity': design_data['boltz_affinity']
        }
        
        # Test first vs last round for each metric that has data
        for metric_name, metric_data in metrics_to_test.items():
            if (len(metric_data) >= 2 and 
                len(metric_data[0]) > 0 and 
                len(metric_data[-1]) > 0):
                
                try:
                    t_stat, p_value = stats.ttest_ind(
                        metric_data[0],
                        metric_data[-1],
                        equal_var=False  # Welch's t-test
                    )
                    
                    significance_annotations[metric_name] = [{
                        'groups': (round_names[0], round_names[-1]),
                        'p_value': p_value
                    }]
                    
                    # Print significance for reference
                    stars = format_significance_label(p_value)
                    print(f"  {metric_name}: {round_names[0]} vs {round_names[-1]} = {stars} (p={p_value:.4f})")
                    
                except Exception as e:
                    print(f"  Could not calculate significance for {metric_name}: {e}")
        
        # Optional: Add pairwise comparisons between consecutive rounds
        for i in range(len(round_names) - 1):
            for metric_name, metric_data in metrics_to_test.items():
                if (len(metric_data) > i+1 and 
                    len(metric_data[i]) > 0 and 
                    len(metric_data[i+1]) > 0):
                    
                    try:
                        t_stat, p_value = stats.ttest_ind(
                            metric_data[i],
                            metric_data[i+1],
                            equal_var=False
                        )
                        
                        # Only add if significant (p < 0.05)
                        if p_value < 0.05:
                            if metric_name not in significance_annotations:
                                significance_annotations[metric_name] = []
                            
                            significance_annotations[metric_name].append({
                                'groups': (round_names[i], round_names[i+1]),
                                'p_value': p_value
                            })
                    except:
                        pass
    
    # Generate all figures with significance annotations
    print("\n" + "="*60)
    print("Generating publication figures with all Boltz metrics...")
    print("="*60)
    
    plot_all_metrics(design_data, args.output_prefix, 
                    metric_significance_annotations=significance_annotations)
    plot_mutation_heatmap(design_data, args.output_prefix)
    
    # Export all data
    export_all_data(design_data, args.output_prefix)
    
    # Print comprehensive statistics
    print_statistics(design_data)
    
    print("\n" + "="*60)
    print("Files generated:")
    print(f"  - {args.output_prefix}_ligandmpnn_confidence.tiff/pdf")
    print(f"  - {args.output_prefix}_seq_rec.tiff/pdf")
    print(f"  - {args.output_prefix}_boltz_confidence.tiff/pdf")
    print(f"  - {args.output_prefix}_boltz_plddt.tiff/pdf")
    print(f"  - {args.output_prefix}_boltz_ptm_iptm.tiff/pdf")
    if any(len(pae) > 0 for pae in design_data['boltz_pocket_pae']):
        print(f"  - {args.output_prefix}_pocket_pae.tiff/pdf")
        print(f"  - {args.output_prefix}_mean_pae.tiff/pdf")
    if any(len(aff) > 0 for aff in design_data['boltz_affinity']):
        print(f"  - {args.output_prefix}_affinity.tiff/pdf")
    if any(len(lig) > 0 for lig in design_data['boltz_ligand_plddt']):
        print(f"  - {args.output_prefix}_chain_plddt.tiff/pdf")
    if any(design_data['boltz_pair_chain_iptm']):
        print(f"  - {args.output_prefix}_pair_chain_iptm.tiff/pdf")
    print(f"  - {args.output_prefix}_composite.tiff/pdf")
    print(f"  - {args.output_prefix}_all_data.csv")
    print("="*60)


if __name__ == "__main__":
    main()