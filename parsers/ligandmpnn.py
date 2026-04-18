#!/usr/bin/env python3
"""
LigandMPNN parser for .fa files with confidence scores
"""

import re
from pathlib import Path


def parse_ligandmpnn_fa(file_path):
    """
    Parse LigandMPNN .fa file format that contains confidence scores in headers
    
    Example header:
    >3mq7_90_48_model_0, id=1, T=0.2, seed=56739, overall_confidence=0.5692, ligand_confidence=0.5692, seq_rec=0.5000
    
    Returns list of dictionaries with sequence data and confidence scores
    """
    sequences = []
    file_path = Path(file_path)
    
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Split by > character to get each entry
        entries = content.split('>')[1:]  # Skip empty first element
        
        for entry in entries:
            lines = entry.strip().split('\n')
            if len(lines) >= 2:
                header = lines[0].strip()
                sequence = ''.join(lines[1:]).replace('\n', '').strip()
                
                # Parse header for confidence scores
                seq_data = {
                    'header': header,
                    'sequence': sequence,
                    'id': None,
                    'overall_confidence': None,
                    'ligand_confidence': None,
                    'seq_rec': None,
                    'temperature': None,
                    'seed': None
                }
                
                # Extract fields using regex
                # overall_confidence
                match = re.search(r'overall_confidence=([0-9.]+)', header)
                if match:
                    seq_data['overall_confidence'] = float(match.group(1))
                
                # ligand_confidence
                match = re.search(r'ligand_confidence=([0-9.]+)', header)
                if match:
                    seq_data['ligand_confidence'] = float(match.group(1))
                
                # seq_rec
                match = re.search(r'seq_rec=([0-9.]+)', header)
                if match:
                    seq_data['seq_rec'] = float(match.group(1))
                
                # id
                match = re.search(r'id=(\d+)', header)
                if match:
                    seq_data['id'] = int(match.group(1))
                
                # temperature
                match = re.search(r'T=([0-9.]+)', header)
                if match:
                    seq_data['temperature'] = float(match.group(1))
                
                # seed
                match = re.search(r'seed=(\d+)', header)
                if match:
                    seq_data['seed'] = int(match.group(1))
                
                sequences.append(seq_data)
    
    except Exception as e:
        print(f"      Error parsing {file_path}: {e}")
    
    return sequences