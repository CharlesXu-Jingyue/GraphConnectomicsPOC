"""
Data validation and basic characterization functions for connectome analysis.
"""

import numpy as np
import pandas as pd
from scipy import stats

def validate_connectivity_matrix(W):
    """
    Validate basic properties of connectivity matrix.
    
    Args:
        W: Connectivity matrix (numpy array)
    
    Returns:
        dict: Validation results
    """
    results = {}
    
    # Basic properties
    results['shape'] = W.shape
    results['is_square'] = W.shape[0] == W.shape[1]
    results['has_negative'] = np.any(W < 0)
    results['has_self_connections'] = np.any(np.diag(W) != 0)
    
    # Sparsity analysis
    nnz = np.count_nonzero(W)
    total_entries = W.size
    results['sparsity'] = 1 - (nnz / total_entries)
    results['density'] = nnz / total_entries
    results['total_connections'] = nnz
    
    # Symmetry assessment
    if results['is_square']:
        symmetry_diff = np.linalg.norm(W - W.T, 'fro') / np.linalg.norm(W, 'fro')
        results['asymmetry_ratio'] = symmetry_diff
        results['is_symmetric'] = symmetry_diff < 1e-10
    
    return results

def compute_degree_statistics(W):
    """
    Compute in-degree and out-degree statistics.
    
    Args:
        W: Connectivity matrix (numpy array)
    
    Returns:
        dict: Degree statistics
    """
    # In-degree: sum of columns (incoming connections)
    in_degree = np.sum(W, axis=0)
    # Out-degree: sum of rows (outgoing connections)  
    out_degree = np.sum(W, axis=1)
    
    results = {
        'in_degree': in_degree,
        'out_degree': out_degree,
        'in_degree_stats': {
            'mean': np.mean(in_degree),
            'std': np.std(in_degree),
            'max': np.max(in_degree),
            'min': np.min(in_degree)
        },
        'out_degree_stats': {
            'mean': np.mean(out_degree),
            'std': np.std(out_degree),
            'max': np.max(out_degree),
            'min': np.min(out_degree)
        }
    }
    
    # Degree correlation
    if len(in_degree) == len(out_degree):
        correlation, p_value = stats.pearsonr(in_degree, out_degree)
        results['degree_correlation'] = correlation
        results['degree_correlation_pvalue'] = p_value
    
    return results

def compute_weight_statistics(W):
    """
    Analyze distribution of connection weights.
    
    Args:
        W: Connectivity matrix (numpy array)
    
    Returns:
        dict: Weight statistics
    """
    # Non-zero weights only
    nonzero_weights = W[W != 0]
    
    if len(nonzero_weights) == 0:
        return {'no_connections': True}
    
    results = {
        'weight_stats': {
            'mean': np.mean(nonzero_weights),
            'std': np.std(nonzero_weights),
            'median': np.median(nonzero_weights),
            'max': np.max(nonzero_weights),
            'min': np.min(nonzero_weights),
            'n_positive': np.sum(nonzero_weights > 0),
            'n_negative': np.sum(nonzero_weights < 0)
        },
        'weights': nonzero_weights
    }
    
    return results