#!/usr/bin/env python3
"""
Pair chain ipTM comparison plotting
"""

import numpy as np
import matplotlib.pyplot as plt

from config import FIG_WIDTH, FIG_HEIGHT
from parsers.boltz import _pair_label_sort


def plot_pair_chain_iptm_comparison(design_data, output_prefix):
    """Compare pair chain ipTM values across rounds with grouped bars."""

    pair_dicts = design_data['boltz_pair_chain_iptm']
    pair_labels = sorted({
        label
        for round_dict in pair_dicts
        for label in round_dict.keys()
    }, key=_pair_label_sort)

    if not pair_labels:
        return

    valid_labels = [
        label for label in pair_labels
        if any(round_dict.get(label) for round_dict in pair_dicts)
    ]
    if not valid_labels:
        return

    round_names = design_data['round_names']
    n_rounds = max(1, len(round_names))
    width = 0.8 / n_rounds
    offsets = (np.arange(n_rounds) - (n_rounds - 1) / 2) * width
    cmap = plt.get_cmap('tab10')
    base_positions = np.arange(len(valid_labels))

    fig, ax = plt.subplots(figsize=(FIG_WIDTH * 1.5, FIG_HEIGHT))

    for round_idx, (round_name, round_dict) in enumerate(zip(round_names, pair_dicts)):
        means = []
        sems = []
        for label in valid_labels:
            values = round_dict.get(label, [])
            if values:
                arr = np.array(values, dtype=float)
                means.append(np.mean(arr))
                sems.append(np.std(arr, ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0)
            else:
                means.append(np.nan)
                sems.append(0)

        valid_idx = [i for i, mean in enumerate(means) if not np.isnan(mean)]
        if not valid_idx:
            continue

        x_positions = base_positions[valid_idx] + offsets[round_idx]
        valid_means = [means[i] for i in valid_idx]
        valid_sems = [sems[i] for i in valid_idx]

        ax.bar(x_positions, valid_means, width=width, label=round_name,
               color=cmap(round_idx % 10), edgecolor='black', linewidth=0.5, alpha=0.93)
        ax.errorbar(x_positions, valid_means, yerr=valid_sems, fmt='none',
                    ecolor='black', capsize=2, capthick=0.5, linewidth=0.5)

        for idx in valid_idx:
            label = valid_labels[idx]
            values = round_dict.get(label, [])
            if not values:
                continue
            jitter = np.random.normal(0, width * 0.12, len(values))
            ax.scatter(base_positions[idx] + offsets[round_idx] + jitter, values,
                       s=4, c='#333333', alpha=0.25, linewidth=0, zorder=2)

    ax.set_xticks(base_positions)
    ax.set_xticklabels(valid_labels, fontsize=6, rotation=45, ha='right')
    ax.set_ylabel('Pair chain ipTM', fontsize=8)
    ax.set_ylim(0, 1)
    ax.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=0.3)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    legend = ax.legend(fontsize=6, frameon=False, loc='upper left',
                       bbox_to_anchor=(1.02, 1))
    if legend:
        legend.set_zorder(5)

    fig.tight_layout(rect=(0, 0, 0.78, 1))
    plt.savefig(f'{output_prefix}_pair_chain_iptm.tiff', dpi=600, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(f'{output_prefix}_pair_chain_iptm.pdf', bbox_inches='tight')
    plt.close(fig)