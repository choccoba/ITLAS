"""Utility functions for ITLAS"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

STAGE_ORDER = ['NL', 'IT', 'IA', 'AR', 'AC']
STAGE_COLORS = {'NL': '#2ecc71', 'IT': '#3498db', 'IA': '#e74c3c', 'AR': '#f39c12', 'AC': '#9b59b6'}

GENE_ID_MAP = {
    'ENSG00000265972': 'MALAT1',
    'ENSG00000117289': 'TXNIP',
    'ENSG00000127528': 'KLF2',
    'ENSG00000198899': 'MT-ATP6',
    'ENSG00000212907': 'MT-ND4L',
    'ENSG00000198938': 'MT-CO3',
}


def summarize_it_findings():
    return """
    ═══════════════════════════════════════════════════════════════
    ITLAS Key Findings: IT-phase Immunopathogenesis
    ═══════════════════════════════════════════════════════════════
    
    1. NK Cell Collapse (Cluster 15)
       - Log2OR = -5.15 in IT (극심한 감소)
       - Recovery in IA phase (Log2OR = +0.28)
    
    2. IT-Exclusive Populations
       - Cluster 21: 1,465 cells (100% IT) - Mito-high/Quiescence
       - Cluster 23: 597 cells (99.8% IT) - Mito-high  
       - Cluster 25: 117 cells (100% IT) - Naive B cell
    
    3. B→Plasma Differentiation Block
       - Cluster 25 (Naive B): IT-exclusive
       - Cluster 22 (Plasma): Depleted across all disease stages
    
    4. Metabolic Stress Signature
       - Mito genes elevated: MT-ATP6, MT-ND4L, MT-CO3
       - KLF2 upregulation: Quiescence marker
    ═══════════════════════════════════════════════════════════════
    """


def plot_stage_comparison(data, value_col, stage_col='Stage', title=''):
    fig, ax = plt.subplots(figsize=(8, 5))
    order = [s for s in STAGE_ORDER if s in data[stage_col].unique()]
    colors = [STAGE_COLORS.get(s, '#333') for s in order]
    sns.boxplot(data=data, x=stage_col, y=value_col, order=order, palette=colors, ax=ax)
    ax.set_title(title)
    plt.tight_layout()
    return fig
