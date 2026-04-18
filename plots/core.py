#!/usr/bin/env python3
"""
Core plotting functions for metric visualization
"""

import numpy as np
import matplotlib.pyplot as plt

from stats.significance import draw_significance_annotations


def plot_metric_with_samples(ax, x_pos, round_names, values, ylabel, color,
                            ylim=None, points=True, significance_annotations=None,
                            sig_p_value_threshold=0.05):
    """Plot metric with multiple samples per round.

    `significance_annotations` accepts a list of dictionaries keyed by `groups`
    (a two-entry tuple of round names) plus either `p_value` for auto-generated
    stars or `text`/`label` to override what is drawn above the bar pair.
    """
    
    # Calculate statistics for each round
    means = []
    sems = []
    valid_rounds = []
    valid_pos = []
    valid_names = []
    n_samples = []
    
    for i, (pos, name, round_values) in enumerate(zip(x_pos, round_names, values)):
        if len(round_values) > 0:
            # Handle potential nested lists (like per-residue arrays) by flattening
            if isinstance(round_values[0], (list, np.ndarray)):
                flat_values = [item for sublist in round_values for item in sublist]
                means.append(np.mean(flat_values))
                if len(flat_values) > 1:
                    sems.append(np.std(flat_values, ddof=1) / np.sqrt(len(flat_values)))
                else:
                    sems.append(0)
                n_samples.append(len(flat_values))
            else:
                means.append(np.mean(round_values))
                if len(round_values) > 1:
                    sems.append(np.std(round_values, ddof=1) / np.sqrt(len(round_values)))
                else:
                    sems.append(0)
                n_samples.append(len(round_values))
            valid_rounds.append(i)
            valid_pos.append(pos)
            valid_names.append(name)
    
    if not valid_rounds:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        return
    
    # Plot bars
    bars = ax.bar(valid_pos, means, width=0.6, color=color,
                  edgecolor='black', linewidth=0.5, alpha=0.9)
    
    # Add error bars
    ax.errorbar(valid_pos, means, yerr=sems, fmt='none',
                capsize=2, capthick=0.5, ecolor='black', elinewidth=0.5)
    
    # Add individual sample points (if not too many)
    if points and max(n_samples) < 100:  # Only show points if not too many
        for i, round_idx in enumerate(valid_rounds):
            round_values = values[round_idx]
            if len(round_values) > 0:
                # For nested lists, flatten for scatter plot
                if isinstance(round_values[0], (list, np.ndarray)):
                    flat_values = [item for sublist in round_values for item in sublist]
                    x_jitter = np.random.normal(valid_pos[i], 0.02, len(flat_values))
                    ax.scatter(x_jitter, flat_values, s=4, c='#AAAAAA',
                              alpha=0.3, linewidth=0, zorder=1)
                else:
                    x_jitter = np.random.normal(valid_pos[i], 0.02, len(round_values))
                    ax.scatter(x_jitter, round_values, s=8, c='#AAAAAA',
                              alpha=0.5, linewidth=0, zorder=1)
    
    # Add sample size annotations
    for i, (pos, mean, sem, n) in enumerate(zip(valid_pos, means, sems, n_samples)):
        y_pos = mean + sem + (ylim[1]*0.02 if ylim else 0.5)
        ax.text(pos, y_pos, f'n={n}', ha='center', va='bottom', fontsize=6)
    
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xticks(valid_pos)
    ax.set_xticklabels(valid_names, fontsize=7, rotation=45, ha='right')
    
    if ylim:
        ax.set_ylim(ylim)

    if significance_annotations:
        draw_significance_annotations(
            ax, valid_names, valid_pos, means, sems,
            significance_annotations, min_p_value=sig_p_value_threshold
        )
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=0.3)
    ax.set_axisbelow(True)