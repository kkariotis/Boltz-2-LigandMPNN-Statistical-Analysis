#!/usr/bin/env python3
"""
Mutation heatmap plotting
"""

import numpy as np
import matplotlib.pyplot as plt

from config import FIG_WIDTH, FIG_HEIGHT


def compute_mutation_heatmap(design_data):
    """Return mutation frequency per position for each round using round 0 as reference."""
    round_sequences = design_data.get('ligandmpnn_sequences') or []
    if len(round_sequences) == 0 or not round_sequences[0]:
        return None, None, None

    reference = None
    for seq in round_sequences[0]:
        if seq:
            reference = seq
            break
    if not reference:
        return None, None, None

    seq_len = len(reference)
    heatmap_rows = []
    rounds_with_data = []

    for round_idx, sequences in enumerate(round_sequences):
        valid_seqs = [s for s in sequences if s and len(s) == seq_len]
        if not valid_seqs:
            heatmap_rows.append(np.full(seq_len, np.nan))
            continue

        mismatch = np.array([
            [1.0 if seq[pos] != reference[pos] else 0.0 for pos in range(seq_len)]
            for seq in valid_seqs
        ])
        heatmap_rows.append(np.mean(mismatch, axis=0))
        rounds_with_data.append(round_idx)

    heatmap_matrix = np.vstack(heatmap_rows)
    return heatmap_matrix, reference, design_data['round_names']


def plot_mutation_heatmap(design_data, output_prefix="figure"):
    """Plot a heatmap of mutation frequency per residue from round 0 reference."""
    heatmap_matrix, reference, round_names = compute_mutation_heatmap(design_data)
    if heatmap_matrix is None:
        return

    n_rounds, seq_len = heatmap_matrix.shape
    width = max(FIG_WIDTH * 1.5, seq_len * 0.08)
    fig, ax = plt.subplots(figsize=(width, FIG_HEIGHT * 1.1))

    im = ax.imshow(heatmap_matrix * 100, cmap='YlOrRd', aspect='auto',
                   vmin=0, vmax=100, origin='lower')

    ax.set_yticks(np.arange(n_rounds))
    ax.set_yticklabels(round_names, fontsize=7)
    xticks = np.linspace(0, seq_len - 1, min(seq_len, 10), dtype=int)
    ax.set_xticks(xticks)
    ax.set_xticklabels((xticks + 1).tolist(), fontsize=6)
    ax.set_xlabel('Residue position')
    ax.set_ylabel('Round')
    ax.set_title('Mutation frequency vs round-0 reference', fontsize=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('% of sequences mutated', fontsize=6)
    cbar.ax.tick_params(labelsize=6)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_mutation_frequency_heatmap.tiff', dpi=600, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(f'{output_prefix}_mutation_frequency_heatmap.pdf', bbox_inches='tight')
    plt.close(fig)