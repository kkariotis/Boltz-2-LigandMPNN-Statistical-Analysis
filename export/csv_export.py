#!/usr/bin/env python3
"""
CSV export and statistics printing functions
"""

import csv
import numpy as np


def export_all_data(design_data, output_prefix="figure"):
    """Export all parsed data to CSV including all Boltz metrics"""
    
    with open(f'{output_prefix}_all_data.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Comprehensive header with all metrics
        writer.writerow([
            'Round', 'Data_Type', 'Sample_ID',
            'LigandMPNN_Confidence', 'LigandMPNN_Ligand_Confidence', 'Sequence_Recovery',
            'Boltz_Confidence', 'Boltz_pLDDT', 'Boltz_PAE', 'Boltz_pTM', 'Boltz_ipTM',
            'Boltz_Ligand_ipTM', 'Boltz_Protein_ipTM',
            'Boltz_Complex_pLDDT', 'Boltz_Complex_pTM', 'Boltz_Complex_ipLDDT', 'Boltz_Complex_PDE', 'Boltz_Complex_iPDE', 'Boltz_Interface_PAE',
            'Boltz_Ligand_pLDDT', 'Boltz_Protein_pLDDT',
            'Boltz_Chain_pTM_Mean', 'Boltz_Chain_pTM_Min', 'Boltz_Chain_pTM_Max', 'Boltz_Pair_ipTM_Mean',
            'Boltz_Mean_PAE_NPZ', 'Boltz_Median_PAE_NPZ', 'Boltz_Std_PAE_NPZ',
            'Boltz_Min_PAE_NPZ', 'Boltz_Max_PAE_NPZ',
            'Boltz_Pocket_PAE', 'Boltz_Pocket_Median_PAE', 'Boltz_Pocket_Std_PAE',
            'Affinity', 'Sequence'
        ])
        
        for round_idx, round_name in enumerate(design_data['round_names']):
            # LigandMPNN entries
            max_mpnn = max(len(design_data['ligandmpnn_confidence'][round_idx]),
                         len(design_data['ligandmpnn_sequences'][round_idx]))
            
            for sample_idx in range(max_mpnn):
                writer.writerow([
                    round_name,
                    'LigandMPNN',
                    sample_idx + 1,
                    f"{design_data['ligandmpnn_confidence'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['ligandmpnn_confidence'][round_idx]) else '',
                    f"{design_data['ligandmpnn_ligand_confidence'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['ligandmpnn_ligand_confidence'][round_idx]) else '',
                    f"{design_data['ligandmpnn_seq_rec'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['ligandmpnn_seq_rec'][round_idx]) else '',
                    # Empty Boltz columns
                    '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
                    design_data['ligandmpnn_sequences'][round_idx][sample_idx] if sample_idx < len(design_data['ligandmpnn_sequences'][round_idx]) else ''
                ])
            
            # Boltz entries
            max_boltz = max(
                len(design_data['boltz_confidence'][round_idx]),
                len(design_data['boltz_plddt'][round_idx]),
                len(design_data['boltz_ptm'][round_idx]),
                len(design_data['boltz_iptm'][round_idx]),
                len(design_data['boltz_pocket_pae'][round_idx]),
                len(design_data['boltz_affinity'][round_idx])
            )
            
            for sample_idx in range(max_boltz):
                row = [
                    round_name,
                    'Boltz',
                    sample_idx + 1,
                    # LigandMPNN columns empty
                    '', '', '',
                    # Boltz confidence metrics
                    f"{design_data['boltz_confidence'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['boltz_confidence'][round_idx]) else '',
                    f"{design_data['boltz_plddt'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_plddt'][round_idx]) else '',
                    f"{design_data['boltz_pae'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_pae'][round_idx]) else '',
                    f"{design_data['boltz_ptm'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['boltz_ptm'][round_idx]) else '',
                    f"{design_data['boltz_iptm'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['boltz_iptm'][round_idx]) else '',
                    f"{design_data['boltz_ligand_iptm'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['boltz_ligand_iptm'][round_idx]) else '',
                    f"{design_data['boltz_protein_iptm'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['boltz_protein_iptm'][round_idx]) else '',
                    f"{design_data['boltz_complex_plddt'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_complex_plddt'][round_idx]) else '',
                    f"{design_data['boltz_complex_ptm'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['boltz_complex_ptm'][round_idx]) else '',
                    f"{design_data['boltz_complex_iplddt'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_complex_iplddt'][round_idx]) else '',
                    f"{design_data['boltz_complex_pde'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_complex_pde'][round_idx]) else '',
                    f"{design_data['boltz_complex_ipde'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_complex_ipde'][round_idx]) else '',
                    f"{design_data['boltz_interface_pae'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_interface_pae'][round_idx]) else '',
                    f"{design_data['boltz_ligand_plddt'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_ligand_plddt'][round_idx]) else '',
                    f"{design_data['boltz_protein_plddt'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_protein_plddt'][round_idx]) else '',
                    f"{design_data['boltz_chains_ptm_mean'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['boltz_chains_ptm_mean'][round_idx]) else '',
                    f"{design_data['boltz_chains_ptm_min'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['boltz_chains_ptm_min'][round_idx]) else '',
                    f"{design_data['boltz_chains_ptm_max'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['boltz_chains_ptm_max'][round_idx]) else '',
                    f"{design_data['boltz_pair_iptm_mean'][round_idx][sample_idx]:.4f}" if sample_idx < len(design_data['boltz_pair_iptm_mean'][round_idx]) else '',
                    # PAE NPZ metrics
                    f"{design_data['boltz_mean_pae_npz'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_mean_pae_npz'][round_idx]) else '',
                    f"{design_data['boltz_median_pae_npz'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_median_pae_npz'][round_idx]) else '',
                    f"{design_data['boltz_std_pae_npz'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_std_pae_npz'][round_idx]) else '',
                    f"{design_data['boltz_min_pae_npz'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_min_pae_npz'][round_idx]) else '',
                    f"{design_data['boltz_max_pae_npz'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_max_pae_npz'][round_idx]) else '',
                    f"{design_data['boltz_pocket_pae'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_pocket_pae'][round_idx]) else '',
                    f"{design_data['boltz_pocket_median_pae'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_pocket_median_pae'][round_idx]) else '',
                    f"{design_data['boltz_pocket_std_pae'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_pocket_std_pae'][round_idx]) else '',
                    # Affinity
                    f"{design_data['boltz_affinity'][round_idx][sample_idx]:.2f}" if sample_idx < len(design_data['boltz_affinity'][round_idx]) else '',
                    # Sequence column empty for Boltz
                    ''
                ]
                writer.writerow(row)
    
    print(f"\nAll data saved to: {output_prefix}_all_data.csv")


def print_statistics(design_data):
    """Print comprehensive statistics for figure captions including all Boltz metrics"""
    
    print("\n" + "="*60)
    print("COMPREHENSIVE STATISTICS FOR FIGURE CAPTION")
    print("="*60)
    
    for i, round_name in enumerate(design_data['round_names']):
        print(f"\n{round_name}:")
        print("-" * 40)
        
        # LigandMPNN stats
        lc = design_data['ligandmpnn_confidence'][i]
        if lc:
            print(f"  LigandMPNN confidence: {np.mean(lc):.4f} ± {np.std(lc, ddof=1):.4f} (n={len(lc)} sequences)")
        
        sr = design_data['ligandmpnn_seq_rec'][i]
        if sr:
            print(f"  Sequence recovery: {np.mean(sr):.4f} ± {np.std(sr, ddof=1):.4f} (n={len(sr)} sequences)")
        
        # Boltz core metrics
        bc = design_data['boltz_confidence'][i]
        if bc:
            print(f"  Boltz confidence: {np.mean(bc):.4f} ± {np.std(bc, ddof=1):.4f} (n={len(bc)} predictions)")
        
        plddt = design_data['boltz_plddt'][i]
        if plddt:
            print(f"  Mean pLDDT: {np.mean(plddt):.2f} ± {np.std(plddt, ddof=1):.2f} (n={len(plddt)})")
        
        ptm = design_data['boltz_ptm'][i]
        if ptm:
            print(f"  pTM: {np.mean(ptm):.4f} ± {np.std(ptm, ddof=1):.4f} (n={len(ptm)})")
        
        iptm = design_data['boltz_iptm'][i]
        if iptm:
            print(f"  ipTM: {np.mean(iptm):.4f} ± {np.std(iptm, ddof=1):.4f} (n={len(iptm)})")

        lig_iptm = design_data['boltz_ligand_iptm'][i]
        if lig_iptm:
            lig_std = np.std(lig_iptm, ddof=1) if len(lig_iptm) > 1 else 0.0
            print(f"  Ligand ipTM: {np.mean(lig_iptm):.4f} ± {lig_std:.4f} (n={len(lig_iptm)})")

        prot_iptm = design_data['boltz_protein_iptm'][i]
        if prot_iptm:
            prot_std = np.std(prot_iptm, ddof=1) if len(prot_iptm) > 1 else 0.0
            print(f"  Protein ipTM: {np.mean(prot_iptm):.4f} ± {prot_std:.4f} (n={len(prot_iptm)})")

        complex_iplddt = design_data['boltz_complex_iplddt'][i]
        if complex_iplddt:
            iplddt_std = np.std(complex_iplddt, ddof=1) if len(complex_iplddt) > 1 else 0.0
            print(f"  Complex ipLDDT: {np.mean(complex_iplddt):.2f} ± {iplddt_std:.2f} (n={len(complex_iplddt)})")

        complex_pde = design_data['boltz_complex_pde'][i]
        if complex_pde:
            pde_std = np.std(complex_pde, ddof=1) if len(complex_pde) > 1 else 0.0
            print(f"  Complex PDE: {np.mean(complex_pde):.2f} ± {pde_std:.2f} (n={len(complex_pde)})")

        complex_ipde = design_data['boltz_complex_ipde'][i]
        if complex_ipde:
            ipde_std = np.std(complex_ipde, ddof=1) if len(complex_ipde) > 1 else 0.0
            print(f"  Complex iPDE: {np.mean(complex_ipde):.2f} ± {ipde_std:.2f} (n={len(complex_ipde)})")

        chain_ptm_mean = design_data['boltz_chains_ptm_mean'][i]
        if chain_ptm_mean:
            ptm_min_vals = design_data['boltz_chains_ptm_min'][i]
            ptm_max_vals = design_data['boltz_chains_ptm_max'][i]
            ptm_min = np.mean(ptm_min_vals) if ptm_min_vals else 0.0
            ptm_max = np.mean(ptm_max_vals) if ptm_max_vals else 0.0
            ptm_std = np.std(chain_ptm_mean, ddof=1) if len(chain_ptm_mean) > 1 else 0.0
            print(f"  Chain pTM mean: {np.mean(chain_ptm_mean):.4f} ± {ptm_std:.4f} (min {ptm_min:.4f}, max {ptm_max:.4f}, n={len(chain_ptm_mean)})")

        pair_iptm = design_data['boltz_pair_iptm_mean'][i]
        if pair_iptm:
            pair_std = np.std(pair_iptm, ddof=1) if len(pair_iptm) > 1 else 0.0
            print(f"  Pairwise ipTM mean: {np.mean(pair_iptm):.4f} ± {pair_std:.4f} (n={len(pair_iptm)})")
        
        pair_chain_data = design_data['boltz_pair_chain_iptm'][i]
        if pair_chain_data:
            flat_pair_values = [val for values in pair_chain_data.values() for val in values]
            if flat_pair_values:
                pair_chain_std = np.std(flat_pair_values, ddof=1) if len(flat_pair_values) > 1 else 0.0
                print(f"  Pair chain ipTM (all pairs): {np.mean(flat_pair_values):.4f} ± {pair_chain_std:.4f} (n={len(flat_pair_values)})")
        
        # PAE metrics
        mean_pae = design_data['boltz_mean_pae_npz'][i]
        if mean_pae:
            print(f"  Mean PAE: {np.mean(mean_pae):.2f} ± {np.std(mean_pae, ddof=1):.2f} Å (n={len(mean_pae)})")
        
        pocket_pae = design_data['boltz_pocket_pae'][i]
        if pocket_pae:
            print(f"  Pocket PAE: {np.mean(pocket_pae):.2f} ± {np.std(pocket_pae, ddof=1):.2f} Å (n={len(pocket_pae)})")
        
        # Affinity
        aff = design_data['boltz_affinity'][i]
        if aff:
            print(f"  Affinity: {np.mean(aff):.2f} ± {np.std(aff, ddof=1):.2f} (n={len(aff)} predictions)")
        
        # Chain-specific metrics
        lig_plddt = design_data['boltz_ligand_plddt'][i]
        if lig_plddt:
            print(f"  Ligand pLDDT: {np.mean(lig_plddt):.2f} ± {np.std(lig_plddt, ddof=1):.2f} (n={len(lig_plddt)})")
        
        prot_plddt = design_data['boltz_protein_plddt'][i]
        if prot_plddt:
            print(f"  Protein pLDDT: {np.mean(prot_plddt):.2f} ± {np.std(prot_plddt, ddof=1):.2f} (n={len(prot_plddt)})")