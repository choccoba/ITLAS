"""IT Signature Score Module"""

import numpy as np
import pandas as pd
from typing import Dict, Union

IT_EXCLUSIVE_MARKERS = {
    'cluster_21': {
        'genes': ['MT-ATP6', 'KLF2', 'MALAT1', 'MT-ND4L', 'MT-CO3'],
        'logfc': [1.57, 1.74, 1.44, 0.80, 0.74],
        'cell_type': 'Mito-high/Quiescence'
    },
    'cluster_23': {
        'genes': ['MT-ATP6', 'KLF2', 'MT-ND4L', 'NACA'],
        'logfc': [1.47, 1.81, 0.85, 1.02],
        'cell_type': 'Mito-high'
    },
    'cluster_25': {
        'genes': ['CD79A', 'MS4A1', 'HLA-DRA', 'CD74'],
        'logfc': [6.05, 5.97, 5.39, 4.20],
        'cell_type': 'IT-exclusive B cell'
    }
}

NK_COLLAPSE_MARKERS = {
    'genes': ['CST7', 'SRGN', 'CXCR4', 'NR4A2', 'DUSP2'],
    'it_log2or': -5.15
}


def calculate_it_signature_score(adata, method='weighted_mean', return_components=False):
    """Calculate IT phase signature score for each cell."""
    scores = {}
    
    for cluster_name, info in IT_EXCLUSIVE_MARKERS.items():
        genes = info['genes']
        weights = np.array(info['logfc'])
        available = [g for g in genes if g in adata.var_names]
        
        if len(available) == 0:
            scores[cluster_name] = np.zeros(adata.n_obs)
            continue
        
        expr = adata[:, available].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()
        
        avail_weights = weights[:len(available)]
        avail_weights = avail_weights / avail_weights.sum()
        scores[cluster_name] = np.average(expr, axis=1, weights=avail_weights)
    
    # NK collapse score
    nk_genes = NK_COLLAPSE_MARKERS['genes']
    available_nk = [g for g in nk_genes if g in adata.var_names]
    if len(available_nk) > 0:
        nk_expr = adata[:, available_nk].X
        if hasattr(nk_expr, 'toarray'):
            nk_expr = nk_expr.toarray()
        scores['nk_collapse'] = -np.mean(nk_expr, axis=1)
    else:
        scores['nk_collapse'] = np.zeros(adata.n_obs)
    
    # Combined score
    weights = {'cluster_21': 0.25, 'cluster_23': 0.25, 'cluster_25': 0.20, 'nk_collapse': 0.30}
    combined = np.zeros(adata.n_obs)
    for comp, w in weights.items():
        if comp in scores:
            s = scores[comp]
            z = (s - np.mean(s)) / (np.std(s) + 1e-8)
            combined += w * z
    scores['IT_signature'] = combined
    
    return scores if return_components else scores['IT_signature']


def validate_against_ground_truth(adata, stage_col='Stage'):
    """Validate IT signature score against known stage labels."""
    scores = calculate_it_signature_score(adata, return_components=True)
    results = []
    
    for stage in adata.obs[stage_col].unique():
        mask = adata.obs[stage_col] == stage
        results.append({
            'Stage': stage,
            'N_cells': mask.sum(),
            'IT_score_mean': np.mean(scores['IT_signature'][mask]),
            'IT_score_std': np.std(scores['IT_signature'][mask])
        })
    
    return pd.DataFrame(results).sort_values('IT_score_mean', ascending=False)
