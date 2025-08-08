# Linkage Methods in Hierarchical Clustering - Explanation
# This explains the different linkage methods and their characteristics

"""
LINKAGE METHODS EXPLANATION
===========================

Linkage methods determine how the distance between clusters is calculated when merging them
in hierarchical clustering. Each method has different characteristics and is suitable for
different types of data and clustering goals.

1. WARD LINKAGE
===============
- **Definition**: Minimizes the within-cluster variance (sum of squared distances)
- **Formula**: d(A,B) = √[(n_A + n_B) * ||c_A - c_B||² / (n_A + n_B)]
- **Characteristics**:
  * Tends to create clusters of similar sizes
  * Sensitive to outliers
  * Works well when clusters are roughly spherical and similar in size
  * Most commonly used method
- **Best for**: Data with roughly equal-sized, spherical clusters
- **Drawback**: Can be sensitive to noise and outliers

2. COMPLETE LINKAGE (Maximum Linkage)
=====================================
- **Definition**: Uses the maximum distance between any two points in different clusters
- **Formula**: d(A,B) = max{d(a,b) : a ∈ A, b ∈ B}
- **Characteristics**:
  * Tends to create compact, spherical clusters
  * Less sensitive to outliers than single linkage
  * Can create clusters of varying sizes
  * More conservative in merging clusters
- **Best for**: Data where you want compact, well-separated clusters
- **Drawback**: Can break large clusters and be sensitive to noise

3. AVERAGE LINKAGE (UPGMA - Unweighted Pair Group Method with Arithmetic Mean)
==============================================================================
- **Definition**: Uses the average distance between all pairs of points in different clusters
- **Formula**: d(A,B) = (1/|A|*|B|) * Σ d(a,b) for all a ∈ A, b ∈ B
- **Characteristics**:
  * Compromise between single and complete linkage
  * Less sensitive to outliers than single linkage
  * Creates clusters of moderate compactness
  * Good balance of cluster size and shape
- **Best for**: General-purpose clustering when you want balanced results
- **Drawback**: Can be computationally intensive for large datasets

4. SINGLE LINKAGE (Minimum Linkage)
===================================
- **Definition**: Uses the minimum distance between any two points in different clusters
- **Formula**: d(A,B) = min{d(a,b) : a ∈ A, b ∈ B}
- **Characteristics**:
  * Tends to create elongated, chain-like clusters
  * Very sensitive to outliers and noise
  * Can create clusters of very different sizes
  * Most aggressive in merging clusters
- **Best for**: Data with elongated or irregular cluster shapes
- **Drawback**: Can create "chaining" effect and be very sensitive to noise

5. CENTROID LINKAGE
===================
- **Definition**: Uses the distance between cluster centroids (mean points)
- **Formula**: d(A,B) = ||c_A - c_B|| where c_A, c_B are centroids
- **Characteristics**:
  * Similar to average linkage but uses centroids
  * Can create clusters of varying sizes
  * Less sensitive to outliers than single linkage
  * Can be unstable in some cases
- **Best for**: Data where centroid-based distance makes sense
- **Drawback**: Can be less stable than other methods

VISUAL COMPARISON
=================

Consider two clusters A and B:

Single Linkage:     Complete Linkage:    Average Linkage:    Ward:
    A ••••              A ••••              A ••••              A ••••
         •••• B               •••• B               •••• B               •••• B
    ↑                    ↑                    ↑                    ↑
  Min dist            Max dist            Avg dist            Variance

PRACTICAL GUIDELINES
====================

1. **Start with Ward**: Good default choice for most datasets
2. **Use Complete**: When you want compact, well-separated clusters
3. **Use Average**: When you want balanced results
4. **Use Single**: Only when you expect elongated clusters
5. **Use Centroid**: When centroid-based distance is meaningful

FOR YOUR COVID-19 DATA
======================

Given that you're analyzing:
- **Angles**: Structural angles (spherical coordinates)
- **Mutations**: Binary/categorical mutation data
- **Distances**: Spatial distances between residues

**Recommended approach**:
1. **Ward linkage**: Best for your data because:
   - Your features are standardized (similar scales)
   - You likely want clusters of similar sizes
   - Structural data often forms roughly spherical clusters

2. **Complete linkage**: Good alternative if you want:
   - More compact, well-separated clusters
   - Less sensitivity to outliers

3. **Average linkage**: Good compromise if Ward doesn't give satisfactory results

**Avoid Single linkage** for your data because:
- Structural data typically doesn't form chain-like clusters
- Very sensitive to noise in mutation data
- Can create misleading elongated clusters
"""

# Example visualization of different linkage methods
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

def demonstrate_linkage_methods():
    """
    Create a simple demonstration of different linkage methods
    """
    # Create sample data with 3 clusters
    np.random.seed(42)
    
    # Cluster 1: centered at (0, 0)
    cluster1 = np.random.normal(0, 1, (20, 2))
    
    # Cluster 2: centered at (4, 4)
    cluster2 = np.random.normal(4, 1, (20, 2))
    
    # Cluster 3: centered at (8, 0)
    cluster3 = np.random.normal(8, 1, (20, 2))
    
    # Combine all data
    data = np.vstack([cluster1, cluster2, cluster3])
    
    # Create labels for the original clusters
    labels = ['Cluster1'] * 20 + ['Cluster2'] * 20 + ['Cluster3'] * 20
    
    # Test different linkage methods
    linkage_methods = ['ward', 'complete', 'average', 'single', 'centroid']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Comparison of Linkage Methods', fontsize=16)
    
    # Plot original data
    ax = axes[0, 0]
    colors = ['red', 'blue', 'green']
    for i, cluster in enumerate([cluster1, cluster2, cluster3]):
        ax.scatter(cluster[:, 0], cluster[:, 1], c=colors[i], alpha=0.7, s=50)
    ax.set_title('Original Data (3 Clusters)')
    ax.set_xlabel('Feature 1')
    ax.set_ylabel('Feature 2')
    ax.grid(True, alpha=0.3)
    
    # Create dendrograms for each method
    for i, method in enumerate(linkage_methods):
        row = (i + 1) // 3
        col = (i + 1) % 3
        
        try:
            # Calculate linkage
            if method == 'ward':
                linkage_matrix = linkage(data, method=method)
            else:
                linkage_matrix = linkage(data, method=method)
            
            # Create dendrogram
            dendrogram(linkage_matrix, 
                      labels=labels,
                      ax=axes[row, col],
                      orientation='top',
                      leaf_rotation=90,
                      leaf_font_size=8)
            
            axes[row, col].set_title(f'{method.capitalize()} Linkage')
            axes[row, col].set_xlabel('Samples')
            axes[row, col].set_ylabel('Distance')
            
        except Exception as e:
            axes[row, col].text(0.5, 0.5, f'Error with {method}\nlinkage method', 
                               ha='center', va='center', transform=axes[row, col].transAxes)
            axes[row, col].set_title(f'{method.capitalize()} Linkage (Error)')
    
    plt.tight_layout()
    plt.show()
    
    return data, labels

# Uncomment to run the demonstration
# data, labels = demonstrate_linkage_methods()

print("Linkage methods explanation completed!")
print("\nKey takeaways for your COVID-19 data:")
print("1. Ward linkage is likely your best choice")
print("2. Complete linkage for more compact clusters")
print("3. Average linkage as a good compromise")
print("4. Avoid single linkage for structural data") 