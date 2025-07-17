"""
Clustering analysis functions for connectome data.
"""

import numpy as np
from sklearn.cluster import SpectralClustering, KMeans, AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def connectivity_based_clustering(W, method='spectral', n_clusters=None, **kwargs):
    """
    Perform clustering based on connectivity patterns.
    
    Args:
        W: Connectivity matrix (numpy array)
        method: Clustering method ('spectral', 'kmeans', 'hierarchical')
        n_clusters: Number of clusters (if None, will estimate)
        **kwargs: Additional parameters for clustering algorithms
    
    Returns:
        dict: Clustering results
    """
    # Prepare connectivity features
    # Use both in-degree and out-degree patterns
    in_degree = np.sum(W, axis=0)
    out_degree = np.sum(W, axis=1)
    
    # Create feature matrix combining connectivity patterns
    # Option 1: Use full connectivity matrix
    features_full = W.copy()
    
    # Option 2: Use degree-based features
    features_degree = np.column_stack([in_degree, out_degree])
    
    # Option 3: Use row-wise connectivity patterns
    features_patterns = W.copy()
    
    # Standardize features
    scaler = StandardScaler()
    features_full_scaled = scaler.fit_transform(features_full)
    features_degree_scaled = scaler.fit_transform(features_degree)
    
    results = {}
    
    # Estimate number of clusters if not provided
    if n_clusters is None:
        n_clusters = estimate_n_clusters(features_full_scaled)
        results['estimated_n_clusters'] = n_clusters
    
    # Perform clustering with different feature sets
    for feature_name, features in [
        ('full_connectivity', features_full_scaled),
        ('degree_patterns', features_degree_scaled),
        ('connectivity_patterns', features_patterns)
    ]:
        
        if method == 'spectral':
            # For spectral clustering, use affinity matrix
            affinity_matrix = np.abs(W) + np.abs(W.T)  # Symmetrize
            clusterer = SpectralClustering(
                n_clusters=n_clusters, 
                affinity='precomputed',
                **kwargs
            )
            labels = clusterer.fit_predict(affinity_matrix)
            
        elif method == 'kmeans':
            clusterer = KMeans(n_clusters=n_clusters, random_state=42, **kwargs)
            labels = clusterer.fit_predict(features)
            
        elif method == 'hierarchical':
            clusterer = AgglomerativeClustering(n_clusters=n_clusters, **kwargs)
            labels = clusterer.fit_predict(features)
        
        # Compute clustering quality metrics
        if len(set(labels)) > 1:  # Need at least 2 clusters for silhouette
            try:
                silhouette = silhouette_score(features, labels)
            except:
                silhouette = np.nan
        else:
            silhouette = np.nan
            
        results[feature_name] = {
            'labels': labels,
            'n_clusters_found': len(set(labels)),
            'silhouette_score': silhouette,
            'cluster_sizes': np.bincount(labels)
        }
    
    return results

def estimate_n_clusters(X, max_clusters=20):
    """
    Estimate optimal number of clusters using silhouette analysis.
    
    Args:
        X: Feature matrix
        max_clusters: Maximum number of clusters to test
    
    Returns:
        int: Estimated optimal number of clusters
    """
    if X.shape[0] < 4:  # Need at least 4 samples
        return 2
        
    max_clusters = min(max_clusters, X.shape[0] // 2)
    silhouette_scores = []
    
    for k in range(2, max_clusters + 1):
        try:
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            labels = kmeans.fit_predict(X)
            score = silhouette_score(X, labels)
            silhouette_scores.append(score)
        except:
            silhouette_scores.append(-1)
    
    if silhouette_scores:
        best_k = np.argmax(silhouette_scores) + 2  # +2 because we start from k=2
        return best_k
    else:
        return 2

def compare_clusterings(pred_labels, true_labels):
    """
    Compare predicted clustering with ground truth labels.
    
    Args:
        pred_labels: Predicted cluster labels
        true_labels: Ground truth labels
    
    Returns:
        dict: Comparison metrics
    """
    # Convert string labels to numeric if needed
    if isinstance(true_labels[0], str):
        unique_labels = list(set(true_labels))
        true_labels_numeric = [unique_labels.index(label) for label in true_labels]
    else:
        true_labels_numeric = true_labels
    
    results = {
        'adjusted_rand_score': adjusted_rand_score(true_labels_numeric, pred_labels),
        'normalized_mutual_info': normalized_mutual_info_score(true_labels_numeric, pred_labels),
        'n_pred_clusters': len(set(pred_labels)),
        'n_true_clusters': len(set(true_labels_numeric))
    }
    
    return results

def analyze_cluster_composition(pred_labels, true_labels):
    """
    Analyze the composition of predicted clusters in terms of true cell types.
    
    Args:
        pred_labels: Predicted cluster labels  
        true_labels: Ground truth cell type labels
    
    Returns:
        dict: Cluster composition analysis
    """
    results = {}
    
    # Create contingency table
    pred_unique = sorted(set(pred_labels))
    true_unique = sorted(set(true_labels))
    
    contingency = np.zeros((len(pred_unique), len(true_unique)))
    
    for i, pred_cluster in enumerate(pred_unique):
        for j, true_type in enumerate(true_unique):
            count = sum((np.array(pred_labels) == pred_cluster) & 
                       (np.array(true_labels) == true_type))
            contingency[i, j] = count
    
    results['contingency_matrix'] = contingency
    results['pred_cluster_labels'] = pred_unique
    results['true_type_labels'] = true_unique
    
    # Analyze purity and recall for each cluster
    cluster_analysis = {}
    for i, pred_cluster in enumerate(pred_unique):
        cluster_mask = np.array(pred_labels) == pred_cluster
        cluster_true_labels = np.array(true_labels)[cluster_mask]
        
        # Most common true type in this cluster
        unique_types, counts = np.unique(cluster_true_labels, return_counts=True)
        most_common_idx = np.argmax(counts)
        most_common_type = unique_types[most_common_idx]
        
        # Purity: fraction of dominant type in cluster
        purity = counts[most_common_idx] / len(cluster_true_labels)
        
        cluster_analysis[pred_cluster] = {
            'size': len(cluster_true_labels),
            'dominant_type': most_common_type,
            'purity': purity,
            'type_distribution': dict(zip(unique_types, counts))
        }
    
    results['cluster_analysis'] = cluster_analysis
    
    return results