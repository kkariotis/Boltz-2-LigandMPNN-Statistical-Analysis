#!/usr/bin/env python3
"""
Individual metric plotting functions
"""

import numpy as np
import matplotlib.pyplot as plt

from config import FIG_WIDTH, FIG_HEIGHT
from plots.core import plot_metric_with_samples
from plots.pair_chain import plot_pair_chain_iptm_comparison


def plot_all_metrics(design_data, output_prefix="figure", metric_significance_annotations=None):
    """Generate all publication figures including all Boltz metrics"""
    
    x_pos = np.arange(len(design_data['round_names']))
    metric_significance_annotations = metric_significance_annotations or {}
    # Expect {'metric_key': [{'groups': ('Round A', 'Round B'), 'p_value': 0.03}, ...], ...}

    def _calc_series_ylim(series, min_margin=0.2):
        flat = [val for values in series for val in values if isinstance(val, (int, float))]
        if not flat:
            return None
        data_min, data_max = min(flat), max(flat)
        delta = data_max - data_min
        if delta == 0:
            margin = max(abs(data_min), abs(data_max), min_margin) * 0.1
        else:
            margin = delta * 0.1
        margin = max(margin, min_margin)
        return (data_min - margin, data_max + margin)
    
    # Figure 1: LigandMPNN Overall Confidence
    fig1, ax1 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    plot_metric_with_samples(
        ax1, x_pos, design_data['round_names'],
        design_data['ligandmpnn_confidence'],
        'LigandMPNN overall confidence',
        '#4C72B0',
        ylim=(0, 1),
        significance_annotations=metric_significance_annotations.get('ligandmpnn_confidence')
    )
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_ligandmpnn_confidence.tiff', dpi=600, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(f'{output_prefix}_ligandmpnn_confidence.pdf', bbox_inches='tight')
    
    # Figure 2: LigandMPNN Sequence Recovery
    fig2, ax2 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    plot_metric_with_samples(
        ax2, x_pos, design_data['round_names'],
        design_data['ligandmpnn_seq_rec'],
        'Sequence recovery',
        '#55A868',
        ylim=(0, 1),
        significance_annotations=metric_significance_annotations.get('ligandmpnn_seq_rec')
    )
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_seq_rec.tiff', dpi=600, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(f'{output_prefix}_seq_rec.pdf', bbox_inches='tight')
    
    # Figure 3: Boltz Confidence Score
    fig3, ax3 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    plot_metric_with_samples(
        ax3, x_pos, design_data['round_names'],
        design_data['boltz_confidence'],
        'Boltz-2 confidence',
        '#C44E52',
        ylim=(0, 1),
        significance_annotations=metric_significance_annotations.get('boltz_confidence')
    )
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_boltz_confidence.tiff', dpi=600, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(f'{output_prefix}_boltz_confidence.pdf', bbox_inches='tight')
    
    # Figure 4: Boltz pLDDT
    if any(len(plddt) > 0 for plddt in design_data['boltz_plddt']):
        fig4, ax4 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
        plot_metric_with_samples(
            ax4, x_pos, design_data['round_names'],
            design_data['boltz_plddt'],
            'Mean pLDDT',
            '#E19C44',
            ylim=(0, 100),
            significance_annotations=metric_significance_annotations.get('boltz_plddt')
        )
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_boltz_plddt.tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        plt.savefig(f'{output_prefix}_boltz_plddt.pdf', bbox_inches='tight')
    
    # Figure 5: pTM and ipTM scores
    if any(len(ptm) > 0 for ptm in design_data['boltz_ptm']) or any(len(iptm) > 0 for iptm in design_data['boltz_iptm']):
        fig5, (ax5a, ax5b) = plt.subplots(1, 2, figsize=(FIG_WIDTH*2, FIG_HEIGHT))
        
        if any(len(ptm) > 0 for ptm in design_data['boltz_ptm']):
            plot_metric_with_samples(
                ax5a, x_pos, design_data['round_names'],
                design_data['boltz_ptm'],
                'pTM score',
                '#4C9AB0',
                ylim=(0, 1),
                significance_annotations=metric_significance_annotations.get('boltz_ptm')
            )
        else:
            ax5a.text(0.5, 0.5, 'No pTM data', ha='center', va='center', transform=ax5a.transAxes)
        
        if any(len(iptm) > 0 for iptm in design_data['boltz_iptm']):
            plot_metric_with_samples(
                ax5b, x_pos, design_data['round_names'],
                design_data['boltz_iptm'],
                'ipTM score',
                '#A14C9A',
                ylim=(0, 1),
                significance_annotations=metric_significance_annotations.get('boltz_iptm')
            )
        else:
            ax5b.text(0.5, 0.5, 'No ipTM data', ha='center', va='center', transform=ax5b.transAxes)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_boltz_ptm_iptm.tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        plt.savefig(f'{output_prefix}_boltz_ptm_iptm.pdf', bbox_inches='tight')
    
    # Figure 6: Pocket PAE
    if any(len(pae) > 0 for pae in design_data['boltz_pocket_pae']):
        fig6, ax6 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
        plot_metric_with_samples(
            ax6, x_pos, design_data['round_names'],
            design_data['boltz_pocket_pae'],
            'Quinone pocket PAE (Å)',
            '#8172B2',
            ylim=(0, 30),
            significance_annotations=metric_significance_annotations.get('boltz_pocket_pae')
        )
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_pocket_pae.tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        plt.savefig(f'{output_prefix}_pocket_pae.pdf', bbox_inches='tight')
    
    # Figure 7: Mean PAE from NPZ
    if any(len(pae) > 0 for pae in design_data['boltz_mean_pae_npz']):
        fig7, ax7 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
        plot_metric_with_samples(
            ax7, x_pos, design_data['round_names'],
            design_data['boltz_mean_pae_npz'],
            'Mean PAE (Å)',
            '#B25A72',
            ylim=(0, 30),
            significance_annotations=metric_significance_annotations.get('boltz_mean_pae_npz')
        )
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_mean_pae.tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        plt.savefig(f'{output_prefix}_mean_pae.pdf', bbox_inches='tight')
    
    # Figure 8: Affinity scores
    if any(len(aff) > 0 for aff in design_data['boltz_affinity']):
        affinity_ylim = _calc_series_ylim(design_data['boltz_affinity'])
        fig8, ax8 = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
        plot_metric_with_samples(
            ax8, x_pos, design_data['round_names'],
            design_data['boltz_affinity'],
            'Binding affinity',
            '#6B5B95',
            ylim=affinity_ylim,
            points=False,
            significance_annotations=metric_significance_annotations.get('boltz_affinity')
        )
        ax8.axhline(0, color='#444444', linewidth=0.6, linestyle=(0, (3, 3)), zorder=0)
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_affinity.tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        plt.savefig(f'{output_prefix}_affinity.pdf', bbox_inches='tight')
    
    # Figure 9: Chain-specific pLDDT (if available)
    if any(len(lig) > 0 for lig in design_data['boltz_ligand_plddt']) or any(len(prot) > 0 for prot in design_data['boltz_protein_plddt']):
        fig9, (ax9a, ax9b) = plt.subplots(1, 2, figsize=(FIG_WIDTH*2, FIG_HEIGHT))
        
        if any(len(lig) > 0 for lig in design_data['boltz_ligand_plddt']):
            plot_metric_with_samples(
                ax9a, x_pos, design_data['round_names'],
                design_data['boltz_ligand_plddt'],
                'Ligand pLDDT',
                '#E19C44',
                ylim=(0, 100),
                significance_annotations=metric_significance_annotations.get('boltz_ligand_plddt')
            )
        else:
            ax9a.text(0.5, 0.5, 'No ligand pLDDT', ha='center', va='center', transform=ax9a.transAxes)
        
        if any(len(prot) > 0 for prot in design_data['boltz_protein_plddt']):
            plot_metric_with_samples(
                ax9b, x_pos, design_data['round_names'],
                design_data['boltz_protein_plddt'],
                'Protein pLDDT',
                '#E19C44',
                ylim=(0, 100),
                significance_annotations=metric_significance_annotations.get('boltz_protein_plddt')
            )
        else:
            ax9b.text(0.5, 0.5, 'No protein pLDDT', ha='center', va='center', transform=ax9b.transAxes)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_chain_plddt.tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        plt.savefig(f'{output_prefix}_chain_plddt.pdf', bbox_inches='tight')

    # Figure 10: Ligand and protein ipTM
    has_ligand_iptm = any(len(vals) > 0 for vals in design_data['boltz_ligand_iptm'])
    has_protein_iptm = any(len(vals) > 0 for vals in design_data['boltz_protein_iptm'])
    if has_ligand_iptm or has_protein_iptm:
        fig_ip, (ax_ip_lig, ax_ip_prot) = plt.subplots(1, 2, figsize=(FIG_WIDTH*2, FIG_HEIGHT))

        if has_ligand_iptm:
            plot_metric_with_samples(
                ax_ip_lig, x_pos, design_data['round_names'],
                design_data['boltz_ligand_iptm'],
                'Ligand ipTM',
                '#4C9AB0',
                ylim=(0, 1),
                significance_annotations=metric_significance_annotations.get('boltz_ligand_iptm')
            )
        else:
            ax_ip_lig.text(0.5, 0.5, 'No ligand ipTM data', ha='center', va='center', transform=ax_ip_lig.transAxes)

        if has_protein_iptm:
            plot_metric_with_samples(
                ax_ip_prot, x_pos, design_data['round_names'],
                design_data['boltz_protein_iptm'],
                'Protein ipTM',
                '#A14C9A',
                ylim=(0, 1),
                significance_annotations=metric_significance_annotations.get('boltz_protein_iptm')
            )
        else:
            ax_ip_prot.text(0.5, 0.5, 'No protein ipTM data', ha='center', va='center', transform=ax_ip_prot.transAxes)

        plt.tight_layout()
        plt.savefig(f'{output_prefix}_chain_iptm.tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        plt.savefig(f'{output_prefix}_chain_iptm.pdf', bbox_inches='tight')

    # Figure 11: Complex-level ipLDDT, PDE, and iPDE
    has_complex_iplddt = any(len(vals) > 0 for vals in design_data['boltz_complex_iplddt'])
    has_complex_pde = any(len(vals) > 0 for vals in design_data['boltz_complex_pde'])
    has_complex_ipde = any(len(vals) > 0 for vals in design_data['boltz_complex_ipde'])
    if has_complex_iplddt or has_complex_pde or has_complex_ipde:
        fig_complex, axes_complex = plt.subplots(1, 3, figsize=(FIG_WIDTH*3, FIG_HEIGHT))
        ax_c_iplddt, ax_c_pde, ax_c_ipde = axes_complex

        if has_complex_iplddt:
            plot_metric_with_samples(
                ax_c_iplddt, x_pos, design_data['round_names'],
                design_data['boltz_complex_iplddt'],
                'Complex ipLDDT',
                '#E19C44',
                ylim=(0, 100),
                significance_annotations=metric_significance_annotations.get('boltz_complex_iplddt')
            )
        else:
            ax_c_iplddt.text(0.5, 0.5, 'No complex ipLDDT', ha='center', va='center', transform=ax_c_iplddt.transAxes)

        if has_complex_pde:
            plot_metric_with_samples(
                ax_c_pde, x_pos, design_data['round_names'],
                design_data['boltz_complex_pde'],
                'Complex PDE',
                '#7F7F7F',
                significance_annotations=metric_significance_annotations.get('boltz_complex_pde')
            )
        else:
            ax_c_pde.text(0.5, 0.5, 'No complex PDE', ha='center', va='center', transform=ax_c_pde.transAxes)

        if has_complex_ipde:
            plot_metric_with_samples(
                ax_c_ipde, x_pos, design_data['round_names'],
                design_data['boltz_complex_ipde'],
                'Complex iPDE',
                '#1F77B4',
                significance_annotations=metric_significance_annotations.get('boltz_complex_ipde')
            )
        else:
            ax_c_ipde.text(0.5, 0.5, 'No complex iPDE', ha='center', va='center', transform=ax_c_ipde.transAxes)

        plt.tight_layout()
        plt.savefig(f'{output_prefix}_complex_iplddt_pde.tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        plt.savefig(f'{output_prefix}_complex_iplddt_pde.pdf', bbox_inches='tight')

    # Figure 12: Chain pTM mean and pairwise ipTM mean
    has_chain_ptm = any(len(vals) > 0 for vals in design_data['boltz_chains_ptm_mean'])
    has_pair_iptm = any(len(vals) > 0 for vals in design_data['boltz_pair_iptm_mean'])
    if has_chain_ptm or has_pair_iptm:
        fig_pair, (ax_pair_ptm, ax_pair_iptm) = plt.subplots(1, 2, figsize=(FIG_WIDTH*2, FIG_HEIGHT))

        if has_chain_ptm:
            plot_metric_with_samples(
                ax_pair_ptm, x_pos, design_data['round_names'],
                design_data['boltz_chains_ptm_mean'],
                'Chain pTM mean',
                '#2CA02C',
                ylim=(0, 1),
                significance_annotations=metric_significance_annotations.get('boltz_chains_ptm_mean')
            )
        else:
            ax_pair_ptm.text(0.5, 0.5, 'No chain pTM mean', ha='center', va='center', transform=ax_pair_ptm.transAxes)

        if has_pair_iptm:
            plot_metric_with_samples(
                ax_pair_iptm, x_pos, design_data['round_names'],
                design_data['boltz_pair_iptm_mean'],
                'Pairwise ipTM mean',
                '#D62728',
                ylim=(0, 1),
                significance_annotations=metric_significance_annotations.get('boltz_pair_iptm_mean')
            )
        else:
            ax_pair_iptm.text(0.5, 0.5, 'No pairwise ipTM', ha='center', va='center', transform=ax_pair_iptm.transAxes)

        plt.tight_layout()
        plt.savefig(f'{output_prefix}_chain_ptm_pair_iptm.tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        plt.savefig(f'{output_prefix}_chain_ptm_pair_iptm.pdf', bbox_inches='tight')
    
    plot_pair_chain_iptm_comparison(design_data, output_prefix)

    # Composite figure with key metrics (4 panels)
    fig_composite, axes = plt.subplots(2, 2, figsize=(7, 6))
    
    metrics = [
        ('ligandmpnn_confidence', design_data['ligandmpnn_confidence'], 'LigandMPNN\nconfidence', '#4C72B0', (0, 1)),
        ('ligandmpnn_seq_rec', design_data['ligandmpnn_seq_rec'], 'Sequence\nrecovery', '#55A868', (0, 1)),
        ('boltz_confidence', design_data['boltz_confidence'], 'Boltz-2\nconfidence', '#C44E52', (0, 1)),
    ]
    
    # Add pocket PAE or affinity as fourth panel
    if any(len(pae) > 0 for pae in design_data['boltz_pocket_pae']):
        metrics.append(('boltz_pocket_pae', design_data['boltz_pocket_pae'], 'Pocket PAE (Å)', '#8172B2', (0, 30)))
    elif any(len(aff) > 0 for aff in design_data['boltz_affinity']):
        metrics.append(('boltz_affinity', design_data['boltz_affinity'], 'Affinity', '#8172B2', None))
    elif any(len(ptm) > 0 for ptm in design_data['boltz_ptm']):
        metrics.append(('boltz_ptm', design_data['boltz_ptm'], 'pTM score', '#4C9AB0', (0, 1)))
    
    for ax, (metric_key, values, ylabel, color, ylim) in zip(axes.flatten(), metrics):
        plot_metric_with_samples(
            ax, x_pos, design_data['round_names'],
            values, ylabel, color, ylim=ylim, points=False,
            significance_annotations=metric_significance_annotations.get(metric_key)
        )
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_composite.tiff', dpi=600, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(f'{output_prefix}_composite.pdf', bbox_inches='tight')
    
    plt.show()