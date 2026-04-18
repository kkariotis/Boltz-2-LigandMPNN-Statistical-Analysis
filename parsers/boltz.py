#!/usr/bin/env python3
"""
Boltz-2 parser for confidence JSON, affinity JSON, and PAE npz files
"""

import json
import numpy as np


def _chain_sort_key(chain_key):
    """Helper for consistent ordering of numeric or string chain identifiers."""
    try:
        return (0, int(chain_key))
    except Exception:
        return (1, str(chain_key))


def _pair_label_sort(label):
    """Sort pair labels by their constituent chain identifiers."""
    parts = label.split('-', 1)
    if len(parts) == 2:
        return (_chain_sort_key(parts[0]), _chain_sort_key(parts[1]))
    return (_chain_sort_key(label),)


def parse_boltz_confidence_json(file_path):
    """Parse Boltz-2 confidence JSON file and extract all metrics"""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        metrics = {
            'confidence_score': None,
            'plddt': None,
            'plddt_per_residue': None,
            'pae': None,
            'pae_matrix': None,
            'ptm': None,
            'iptm': None,
            'ligand_iptm': None,
            'protein_iptm': None,
            'complex_plddt': None,
            'complex_ptm': None,
            'complex_iplddt': None,
            'complex_pde': None,
            'complex_ipde': None,
            'interface_pae': None,
            'ligand_plddt': None,
            'protein_plddt': None,
            'chains_ptm_mean': None,
            'chains_ptm_min': None,
            'chains_ptm_max': None,
            'pair_iptm_mean': None,
            'pair_chain_iptm': None
        }
        
        # Extract overall confidence
        if 'confidence_score' in data:
            metrics['confidence_score'] = data['confidence_score']
        
        # Extract pLDDT if available
        if 'plddt' in data:
            if isinstance(data['plddt'], list):
                metrics['plddt'] = np.mean(data['plddt'])
                metrics['plddt_per_residue'] = data['plddt']
            else:
                metrics['plddt'] = data['plddt']
        
        # Extract PAE info
        if 'pae' in data:
            if isinstance(data['pae'], list):
                pae_matrix = np.array(data['pae'])
                metrics['pae'] = np.mean(pae_matrix)
                metrics['pae_matrix'] = pae_matrix.tolist()
            else:
                metrics['pae'] = data['pae']
        
        # Extract pTM and ipTM
        if 'ptm' in data:
            metrics['ptm'] = data['ptm']
        if 'iptm' in data:
            metrics['iptm'] = data['iptm']
        if 'ligand_iptm' in data:
            metrics['ligand_iptm'] = data['ligand_iptm']
        if 'protein_iptm' in data:
            metrics['protein_iptm'] = data['protein_iptm']
        
        # Extract complex-level metrics
        if 'complex_plddt' in data:
            metrics['complex_plddt'] = data['complex_plddt']
        if 'complex_ptm' in data:
            metrics['complex_ptm'] = data['complex_ptm']
        if 'complex_iplddt' in data:
            metrics['complex_iplddt'] = data['complex_iplddt']
        if 'complex_pde' in data:
            metrics['complex_pde'] = data['complex_pde']
        if 'complex_ipde' in data:
            metrics['complex_ipde'] = data['complex_ipde']
        if 'interface_pae' in data:
            metrics['interface_pae'] = data['interface_pae']
        
        # Extract chain-specific metrics if available
        if 'chain_scores' in data:
            chain_data = data['chain_scores']
            if 'ligand' in chain_data and 'plddt' in chain_data['ligand']:
                metrics['ligand_plddt'] = np.mean(chain_data['ligand']['plddt']) if isinstance(chain_data['ligand']['plddt'], list) else chain_data['ligand']['plddt']
            if 'protein' in chain_data and 'plddt' in chain_data['protein']:
                metrics['protein_plddt'] = np.mean(chain_data['protein']['plddt']) if isinstance(chain_data['protein']['plddt'], list) else chain_data['protein']['plddt']

        if 'chains_ptm' in data and isinstance(data['chains_ptm'], dict):
            chain_values = [float(val) for val in data['chains_ptm'].values() if isinstance(val, (int, float))]
            if chain_values:
                metrics['chains_ptm_mean'] = np.mean(chain_values)
                metrics['chains_ptm_min'] = np.min(chain_values)
                metrics['chains_ptm_max'] = np.max(chain_values)

        if 'pair_chains_iptm' in data and isinstance(data['pair_chains_iptm'], dict):
            pair_map = {}
            pair_values = []
            for row_key in sorted(data['pair_chains_iptm'].keys(), key=_chain_sort_key):
                row_dict = data['pair_chains_iptm'][row_key]
                if not isinstance(row_dict, dict):
                    continue

                for col_key in sorted(row_dict.keys(), key=_chain_sort_key):
                    if _chain_sort_key(row_key) > _chain_sort_key(col_key):
                        continue

                    val_ij = row_dict.get(col_key)
                    if not isinstance(val_ij, (int, float)):
                        continue

                    collected = [float(val_ij)]
                    if row_key != col_key:
                        other_val = data['pair_chains_iptm'].get(col_key, {}).get(row_key)
                        if isinstance(other_val, (int, float)):
                            collected.append(float(other_val))

                    label_keys = sorted([row_key, col_key], key=_chain_sort_key)
                    label = f"{label_keys[0]}-{label_keys[1]}"
                    pair_value = float(np.mean(collected))
                    pair_map[label] = pair_value
                    pair_values.append(pair_value)

            if pair_map:
                metrics['pair_chain_iptm'] = pair_map
            if pair_values:
                metrics['pair_iptm_mean'] = np.mean(pair_values)
        
        return metrics
    
    except Exception as e:
        print(f"      Error parsing Boltz confidence JSON: {e}")
        return None


def parse_boltz_affinity_json(file_path):
    """Parse Boltz-2 affinity JSON file"""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        affinities = []
        def _append_affinity_values(value):
            if isinstance(value, (int, float)):
                affinities.append(value)
            elif isinstance(value, (list, tuple)):
                affinities.extend([v for v in value if isinstance(v, (int, float))])
        
        # Extract affinity scores (structure may vary)
        if isinstance(data, dict):
            for key, value in data.items():
                key_lower = key.lower()
                if key_lower != 'affinity' and 'affinity' in key_lower and 'probability' not in key_lower:
                    _append_affinity_values(value)
            if 'affinity' in data:
                _append_affinity_values(data['affinity'])
            if 'scores' in data:
                _append_affinity_values(data['scores'])
            if 'binding_affinity' in data:
                _append_affinity_values(data['binding_affinity'])
        elif isinstance(data, list):
            _append_affinity_values(data)
        
        return np.mean(affinities) if affinities else None
    
    except Exception as e:
        print(f"      Error parsing Boltz affinity JSON: {e}")
        return None


def parse_boltz_pae_npz(file_path, binding_residues=None):
    """Parse Boltz-2 PAE npz file and extract pocket-specific PAE"""
    try:
        data = np.load(file_path)
        
        # Usually PAE is stored under 'pae' or 'arr_0'
        if 'pae' in data:
            pae_matrix = data['pae']
        else:
            pae_matrix = data['arr_0']
        
        metrics = {
            'mean_pae': np.mean(pae_matrix),
            'median_pae': np.median(pae_matrix),
            'std_pae': np.std(pae_matrix),
            'min_pae': np.min(pae_matrix),
            'max_pae': np.max(pae_matrix),
            'pocket_pae': None,
            'pocket_median_pae': None,
            'pocket_std_pae': None
        }
        
        # Calculate pocket-specific PAE if binding residues provided
        if binding_residues and len(pae_matrix) > 0:
            # Convert 1-indexed residues to 0-indexed
            pocket_indices = [i-1 for i in binding_residues if i-1 < pae_matrix.shape[0]]
            if pocket_indices:
                # Extract submatrix for pocket residues
                pocket_pae = pae_matrix[np.ix_(pocket_indices, pocket_indices)]
                metrics['pocket_pae'] = np.mean(pocket_pae)
                metrics['pocket_median_pae'] = np.median(pocket_pae)
                metrics['pocket_std_pae'] = np.std(pocket_pae)
        
        return metrics
    
    except Exception as e:
        print(f"      Error parsing Boltz PAE npz: {e}")
        return None