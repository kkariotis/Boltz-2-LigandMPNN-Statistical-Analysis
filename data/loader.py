#!/usr/bin/env python3
"""
Main data loader for quinone-binding cytochrome design data
"""

import numpy as np
from pathlib import Path
from collections import defaultdict

from parsers.ligandmpnn import parse_ligandmpnn_fa
from parsers.boltz import (
    parse_boltz_confidence_json,
    parse_boltz_affinity_json,
    parse_boltz_pae_npz
)


def load_quinone_design_data(base_folder, binding_residues=None):
    """
    Load data from round_x folders with recursive file search
    Finds files regardless of nesting depth
    Includes all Boltz confidence metrics
    """
    
    print(f"\n{'='*60}")
    print(f"Loading quinone-binding cytochrome design data from:")
    print(f"{base_folder}")
    print(f"{'='*60}")
    
    base = Path(base_folder)
    
    # Find all round folders
    round_folders = sorted([f for f in base.glob("round_*") if f.is_dir()])
    print(f"\nFound {len(round_folders)} design rounds")
    
    # Initialize data structures with all metrics
    design_data = {
        'round_names': [],
        # LigandMPNN metrics
        'ligandmpnn_confidence': [],
        'ligandmpnn_ligand_confidence': [],
        'ligandmpnn_seq_rec': [],
        'ligandmpnn_sequences': [],
        
        # Boltz confidence metrics (from JSON)
        'boltz_confidence': [],
        'boltz_plddt': [],
        'boltz_plddt_per_residue': [],
        'boltz_pae': [],
        'boltz_pae_matrix': [],
        'boltz_ptm': [],
        'boltz_iptm': [],
        'boltz_complex_plddt': [],
        'boltz_complex_ptm': [],
        'boltz_interface_pae': [],
        'boltz_ligand_iptm': [],
        'boltz_protein_iptm': [],
        'boltz_complex_iplddt': [],
        'boltz_complex_pde': [],
        'boltz_complex_ipde': [],
        'boltz_chains_ptm_mean': [],
        'boltz_chains_ptm_min': [],
        'boltz_chains_ptm_max': [],
        'boltz_pair_iptm_mean': [],
        'boltz_ligand_plddt': [],
        'boltz_protein_plddt': [],
        'boltz_pair_chain_iptm': [],
        
        # Boltz PAE metrics (from NPZ)
        'boltz_mean_pae_npz': [],
        'boltz_median_pae_npz': [],
        'boltz_std_pae_npz': [],
        'boltz_min_pae_npz': [],
        'boltz_max_pae_npz': [],
        'boltz_pocket_pae': [],
        'boltz_pocket_median_pae': [],
        'boltz_pocket_std_pae': [],
        
        # Affinity metrics
        'boltz_affinity': [],
        
        # Sample counts
        'sample_counts': []
    }
    
    for round_idx, round_folder in enumerate(round_folders):
        round_name = round_folder.name.replace('_', ' ').title()
        print(f"\n  {'='*40}")
        print(f"  Processing {round_name}:")
        print(f"  {'='*40}")
        
        # --- Process LigandMPNN .fa files (recursive search) ---
        mpnn_base = round_folder / "ligandmpnn"
        round_mpnn_confidence = []
        round_mpnn_ligand_conf = []
        round_mpnn_seq_rec = []
        round_mpnn_sequences = []
        
        if mpnn_base.exists():
            # Recursively find all .fa files
            fa_files = list(mpnn_base.rglob("*.fa"))
            print(f"    Found {len(fa_files)} LigandMPNN .fa files")
            
            for fa_file in fa_files:
                sequences = parse_ligandmpnn_fa(fa_file)
                
                for seq in sequences:
                    if seq['overall_confidence'] is not None:
                        round_mpnn_confidence.append(seq['overall_confidence'])
                    if seq['ligand_confidence'] is not None:
                        round_mpnn_ligand_conf.append(seq['ligand_confidence'])
                    if seq['seq_rec'] is not None:
                        round_mpnn_seq_rec.append(seq['seq_rec'])
                    if seq['sequence']:
                        round_mpnn_sequences.append(seq['sequence'])
            
            print(f"      Extracted {len(round_mpnn_confidence)} sequences with confidence scores")
            if round_mpnn_confidence:
                print(f"        Overall confidence: {np.mean(round_mpnn_confidence):.4f} ± {np.std(round_mpnn_confidence, ddof=1):.4f}")
                print(f"        Sequence recovery: {np.mean(round_mpnn_seq_rec):.4f} ± {np.std(round_mpnn_seq_rec, ddof=1):.4f}")
        else:
            print(f"    No LigandMPNN folder found at {mpnn_base}")
        
        design_data['ligandmpnn_confidence'].append(round_mpnn_confidence)
        design_data['ligandmpnn_ligand_confidence'].append(round_mpnn_ligand_conf)
        design_data['ligandmpnn_seq_rec'].append(round_mpnn_seq_rec)
        design_data['ligandmpnn_sequences'].append(round_mpnn_sequences)
        
        # --- Process Boltz-2 predictions (recursive search) ---
        boltz_base = round_folder / "boltz"
        
        # Initialize round-specific lists for all Boltz metrics
        round_boltz_confidence = []
        round_boltz_plddt = []
        round_boltz_plddt_per_residue = []
        round_boltz_pae = []
        round_boltz_pae_matrix = []
        round_boltz_ptm = []
        round_boltz_iptm = []
        round_boltz_complex_plddt = []
        round_boltz_complex_ptm = []
        round_boltz_interface_pae = []
        round_boltz_ligand_iptm = []
        round_boltz_protein_iptm = []
        round_boltz_complex_iplddt = []
        round_boltz_complex_pde = []
        round_boltz_complex_ipde = []
        round_boltz_chains_ptm_mean = []
        round_boltz_chains_ptm_min = []
        round_boltz_chains_ptm_max = []
        round_boltz_pair_iptm_mean = []
        round_boltz_ligand_plddt = []
        round_boltz_protein_plddt = []
        round_boltz_pair_chain_iptm = defaultdict(list)
        
        round_boltz_mean_pae_npz = []
        round_boltz_median_pae_npz = []
        round_boltz_std_pae_npz = []
        round_boltz_min_pae_npz = []
        round_boltz_max_pae_npz = []
        round_boltz_pocket_pae = []
        round_boltz_pocket_median_pae = []
        round_boltz_pocket_std_pae = []
        
        round_boltz_affinity = []
        
        if boltz_base.exists():
            # Recursively find all confidence JSON files
            conf_files = list(boltz_base.rglob("confidence_*.json"))
            affinity_files = list(boltz_base.rglob("affinity_*.json"))
            pae_files = list(boltz_base.rglob("pae_*.npz"))
            
            print(f"    Found {len(conf_files)} Boltz confidence files")
            print(f"    Found {len(affinity_files)} Boltz affinity files")
            print(f"    Found {len(pae_files)} Boltz PAE files")
            
            # Process confidence files - extract ALL metrics
            for conf_file in conf_files:
                metrics = parse_boltz_confidence_json(conf_file)
                if metrics:
                    if metrics['confidence_score'] is not None:
                        round_boltz_confidence.append(metrics['confidence_score'])
                    if metrics['plddt'] is not None:
                        round_boltz_plddt.append(metrics['plddt'])
                    if metrics['plddt_per_residue'] is not None:
                        round_boltz_plddt_per_residue.append(metrics['plddt_per_residue'])
                    if metrics['pae'] is not None:
                        round_boltz_pae.append(metrics['pae'])
                    if metrics['pae_matrix'] is not None:
                        round_boltz_pae_matrix.append(metrics['pae_matrix'])
                    if metrics['ptm'] is not None:
                        round_boltz_ptm.append(metrics['ptm'])
                    if metrics['iptm'] is not None:
                        round_boltz_iptm.append(metrics['iptm'])
                    if metrics['complex_plddt'] is not None:
                        round_boltz_complex_plddt.append(metrics['complex_plddt'])
                    if metrics['complex_ptm'] is not None:
                        round_boltz_complex_ptm.append(metrics['complex_ptm'])
                    if metrics['complex_iplddt'] is not None:
                        round_boltz_complex_iplddt.append(metrics['complex_iplddt'])
                    if metrics['complex_pde'] is not None:
                        round_boltz_complex_pde.append(metrics['complex_pde'])
                    if metrics['complex_ipde'] is not None:
                        round_boltz_complex_ipde.append(metrics['complex_ipde'])
                    if metrics['interface_pae'] is not None:
                        round_boltz_interface_pae.append(metrics['interface_pae'])
                    if metrics['ligand_iptm'] is not None:
                        round_boltz_ligand_iptm.append(metrics['ligand_iptm'])
                    if metrics['protein_iptm'] is not None:
                        round_boltz_protein_iptm.append(metrics['protein_iptm'])
                    if metrics['chains_ptm_mean'] is not None:
                        round_boltz_chains_ptm_mean.append(metrics['chains_ptm_mean'])
                    if metrics['chains_ptm_min'] is not None:
                        round_boltz_chains_ptm_min.append(metrics['chains_ptm_min'])
                    if metrics['chains_ptm_max'] is not None:
                        round_boltz_chains_ptm_max.append(metrics['chains_ptm_max'])
                    if metrics['pair_iptm_mean'] is not None:
                        round_boltz_pair_iptm_mean.append(metrics['pair_iptm_mean'])
                    if metrics['pair_chain_iptm']:
                        for label, value in metrics['pair_chain_iptm'].items():
                            round_boltz_pair_chain_iptm[label].append(value)
                    if metrics['ligand_plddt'] is not None:
                        round_boltz_ligand_plddt.append(metrics['ligand_plddt'])
                    if metrics['protein_plddt'] is not None:
                        round_boltz_protein_plddt.append(metrics['protein_plddt'])
            
            # Process affinity files
            for aff_file in affinity_files:
                affinity = parse_boltz_affinity_json(aff_file)
                if affinity is not None:
                    round_boltz_affinity.append(affinity)
            
            # Process PAE files for pocket-specific metrics
            for pae_file in pae_files:
                pae_metrics = parse_boltz_pae_npz(pae_file, binding_residues)
                if pae_metrics:
                    if pae_metrics['mean_pae'] is not None:
                        round_boltz_mean_pae_npz.append(pae_metrics['mean_pae'])
                    if pae_metrics['median_pae'] is not None:
                        round_boltz_median_pae_npz.append(pae_metrics['median_pae'])
                    if pae_metrics['std_pae'] is not None:
                        round_boltz_std_pae_npz.append(pae_metrics['std_pae'])
                    if pae_metrics['min_pae'] is not None:
                        round_boltz_min_pae_npz.append(pae_metrics['min_pae'])
                    if pae_metrics['max_pae'] is not None:
                        round_boltz_max_pae_npz.append(pae_metrics['max_pae'])
                    if pae_metrics['pocket_pae'] is not None:
                        round_boltz_pocket_pae.append(pae_metrics['pocket_pae'])
                    if pae_metrics['pocket_median_pae'] is not None:
                        round_boltz_pocket_median_pae.append(pae_metrics['pocket_median_pae'])
                    if pae_metrics['pocket_std_pae'] is not None:
                        round_boltz_pocket_std_pae.append(pae_metrics['pocket_std_pae'])
            
            # Print comprehensive summary
            print(f"      Boltz confidence metrics extracted:")
            if round_boltz_confidence:
                print(f"        - Confidence score: {np.mean(round_boltz_confidence):.4f} ± {np.std(round_boltz_confidence, ddof=1):.4f} (n={len(round_boltz_confidence)})")
            if round_boltz_plddt:
                print(f"        - Mean pLDDT: {np.mean(round_boltz_plddt):.2f} ± {np.std(round_boltz_plddt, ddof=1):.2f} (n={len(round_boltz_plddt)})")
            if round_boltz_ptm:
                print(f"        - pTM: {np.mean(round_boltz_ptm):.4f} ± {np.std(round_boltz_ptm, ddof=1):.4f} (n={len(round_boltz_ptm)})")
            if round_boltz_iptm:
                print(f"        - ipTM: {np.mean(round_boltz_iptm):.4f} ± {np.std(round_boltz_iptm, ddof=1):.4f} (n={len(round_boltz_iptm)})")
            if round_boltz_ligand_iptm:
                ligand_mean = np.mean(round_boltz_ligand_iptm)
                ligand_std = np.std(round_boltz_ligand_iptm, ddof=1) if len(round_boltz_ligand_iptm) > 1 else 0.0
                print(f"        - Ligand ipTM: {ligand_mean:.4f} ± {ligand_std:.4f} (n={len(round_boltz_ligand_iptm)})")
            if round_boltz_protein_iptm:
                protein_mean = np.mean(round_boltz_protein_iptm)
                protein_std = np.std(round_boltz_protein_iptm, ddof=1) if len(round_boltz_protein_iptm) > 1 else 0.0
                print(f"        - Protein ipTM: {protein_mean:.4f} ± {protein_std:.4f} (n={len(round_boltz_protein_iptm)})")
            if round_boltz_complex_iplddt:
                print(f"        - Complex ipLDDT: {np.mean(round_boltz_complex_iplddt):.2f} ± {np.std(round_boltz_complex_iplddt, ddof=1):.2f} (n={len(round_boltz_complex_iplddt)})")
            if round_boltz_complex_pde:
                print(f"        - Complex PDE: {np.mean(round_boltz_complex_pde):.2f} ± {np.std(round_boltz_complex_pde, ddof=1):.2f} (n={len(round_boltz_complex_pde)})")
            if round_boltz_complex_ipde:
                print(f"        - Complex iPDE: {np.mean(round_boltz_complex_ipde):.2f} ± {np.std(round_boltz_complex_ipde, ddof=1):.2f} (n={len(round_boltz_complex_ipde)})")
            if round_boltz_pocket_pae:
                print(f"        - Pocket PAE: {np.mean(round_boltz_pocket_pae):.2f} ± {np.std(round_boltz_pocket_pae, ddof=1):.2f} Å (n={len(round_boltz_pocket_pae)})")
            if round_boltz_affinity:
                print(f"        - Affinity: {np.mean(round_boltz_affinity):.2f} ± {np.std(round_boltz_affinity, ddof=1):.2f} (n={len(round_boltz_affinity)})")
            if round_boltz_chains_ptm_mean:
                ptm_mean = np.mean(round_boltz_chains_ptm_mean)
                ptm_std = np.std(round_boltz_chains_ptm_mean, ddof=1) if len(round_boltz_chains_ptm_mean) > 1 else 0.0
                ptm_min = np.mean(round_boltz_chains_ptm_min) if round_boltz_chains_ptm_min else 0.0
                ptm_max = np.mean(round_boltz_chains_ptm_max) if round_boltz_chains_ptm_max else 0.0
                print(f"        - Chain pTM mean: {ptm_mean:.4f} ± {ptm_std:.4f} (min {ptm_min:.4f}, max {ptm_max:.4f}, n={len(round_boltz_chains_ptm_mean)})")
            if round_boltz_pair_iptm_mean:
                pair_mean = np.mean(round_boltz_pair_iptm_mean)
                pair_std = np.std(round_boltz_pair_iptm_mean, ddof=1) if len(round_boltz_pair_iptm_mean) > 1 else 0.0
                print(f"        - Pairwise ipTM mean: {pair_mean:.4f} ± {pair_std:.4f} (n={len(round_boltz_pair_iptm_mean)})")
        
        else:
            print(f"    No Boltz folder found at {boltz_base}")
        
        # Store all Boltz metrics
        design_data['boltz_confidence'].append(round_boltz_confidence)
        design_data['boltz_plddt'].append(round_boltz_plddt)
        design_data['boltz_plddt_per_residue'].append(round_boltz_plddt_per_residue)
        design_data['boltz_pae'].append(round_boltz_pae)
        design_data['boltz_pae_matrix'].append(round_boltz_pae_matrix)
        design_data['boltz_ptm'].append(round_boltz_ptm)
        design_data['boltz_iptm'].append(round_boltz_iptm)
        design_data['boltz_complex_plddt'].append(round_boltz_complex_plddt)
        design_data['boltz_complex_ptm'].append(round_boltz_complex_ptm)
        design_data['boltz_interface_pae'].append(round_boltz_interface_pae)
        design_data['boltz_ligand_iptm'].append(round_boltz_ligand_iptm)
        design_data['boltz_protein_iptm'].append(round_boltz_protein_iptm)
        design_data['boltz_complex_iplddt'].append(round_boltz_complex_iplddt)
        design_data['boltz_complex_pde'].append(round_boltz_complex_pde)
        design_data['boltz_complex_ipde'].append(round_boltz_complex_ipde)
        design_data['boltz_chains_ptm_mean'].append(round_boltz_chains_ptm_mean)
        design_data['boltz_chains_ptm_min'].append(round_boltz_chains_ptm_min)
        design_data['boltz_chains_ptm_max'].append(round_boltz_chains_ptm_max)
        design_data['boltz_pair_iptm_mean'].append(round_boltz_pair_iptm_mean)
        design_data['boltz_pair_chain_iptm'].append(dict(round_boltz_pair_chain_iptm))
        design_data['boltz_ligand_plddt'].append(round_boltz_ligand_plddt)
        design_data['boltz_protein_plddt'].append(round_boltz_protein_plddt)
        
        design_data['boltz_mean_pae_npz'].append(round_boltz_mean_pae_npz)
        design_data['boltz_median_pae_npz'].append(round_boltz_median_pae_npz)
        design_data['boltz_std_pae_npz'].append(round_boltz_std_pae_npz)
        design_data['boltz_min_pae_npz'].append(round_boltz_min_pae_npz)
        design_data['boltz_max_pae_npz'].append(round_boltz_max_pae_npz)
        design_data['boltz_pocket_pae'].append(round_boltz_pocket_pae)
        design_data['boltz_pocket_median_pae'].append(round_boltz_pocket_median_pae)
        design_data['boltz_pocket_std_pae'].append(round_boltz_pocket_std_pae)
        
        design_data['boltz_affinity'].append(round_boltz_affinity)
        
        # Store sample counts
        design_data['sample_counts'].append({
            'ligandmpnn': len(round_mpnn_confidence),
            'boltz_confidence': len(round_boltz_confidence),
            'boltz_plddt': len(round_boltz_plddt),
            'boltz_ptm': len(round_boltz_ptm),
            'boltz_iptm': len(round_boltz_iptm),
            'boltz_ligand_iptm': len(round_boltz_ligand_iptm),
            'boltz_complex_iplddt': len(round_boltz_complex_iplddt),
            'boltz_pair_iptm_mean': len(round_boltz_pair_iptm_mean),
            'boltz_pair_chain_iptm': sum(len(vals) for vals in round_boltz_pair_chain_iptm.values()),
            'boltz_pocket_pae': len(round_boltz_pocket_pae),
            'boltz_affinity': len(round_boltz_affinity)
        })
        
        design_data['round_names'].append(round_name)
    
    return design_data