# Comprehensive Clustering Dendrogram Analysis - Notebook Version
# Add this code to your existing notebook to create dendrograms with all three data types

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cophenet
from scipy.spatial.distance import pdist, squareform
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 1. LOAD AND COMBINE ALL DATASETS
# ============================================================================

print("=" * 70)
print("COMPREHENSIVE CLUSTERING DENDROGRAM ANALYSIS")
print("=" * 70)

print("\n1. Loading datasets...")

# Load sidechain angles data
print("   Loading sidechain angles data...")
sidechain_df = pd.read_csv("sidechain_angles.csv", index_col=0)

# Extract variant names from sample names
sidechain_df['variant'] = [sample.split('-S')[0] if '-S' in sample else sample for sample in sidechain_df.index]
sidechain_df = sidechain_df.reset_index().rename(columns={'index': 'sequence'})

# Handle missing values in angles
angle_cols = [col for col in sidechain_df.columns if any(col.endswith(suffix) for suffix in ['_Theta', '_Phi', '_Chi1'])]
sidechain_df[angle_cols] = sidechain_df[angle_cols].replace('N/A', np.nan).astype(float)
sidechain_df[angle_cols] = sidechain_df[angle_cols].fillna(sidechain_df[angle_cols].mean())

print(f"   Sidechain data shape: {sidechain_df.shape}")
print(f"   Angle features: {len(angle_cols)}")

# Load mutation data
print("   Loading mutation data...")
mutation_df = pd.read_csv("mutation/mutation-all.csv")
print(f"   Mutation data shape: {mutation_df.shape}")

# Load distance data
print("   Loading distance data...")
try:
    distance_df = pd.read_csv("tmalign/residues_full.csv")
    print(f"   Distance data shape: {distance_df.shape}")
except FileNotFoundError:
    print("   Distance data not found, proceeding without it")
    distance_df = None

# Combine datasets
print("\n2. Combining datasets...")
combined_df = sidechain_df.copy()

# Merge with mutation data
combined_df = combined_df.merge(mutation_df, on=['sequence', 'variant'], how='left')
print(f"   After merging mutation data: {combined_df.shape}")

# Merge with distance data if available
if distance_df is not None:
    combined_df = combined_df.merge(distance_df, on='sequence', how='left')
    print(f"   After merging distance data: {combined_df.shape}")

# Remove rows with missing data
combined_df = combined_df.dropna()
print(f"   Final combined dataset shape: {combined_df.shape}")

# ============================================================================
# 2. EXTRACT FEATURES FOR CLUSTERING
# ============================================================================

print("\n3. Extracting features for clustering...")

# Define feature groups
feature_groups = {}

# Angle features
feature_groups['angles'] = angle_cols
print(f"   Angle features: {len(angle_cols)}")

# Mutation features (exclude sequence, variant, and 'Against' columns)
mutation_cols = [col for col in combined_df.columns 
                if col not in ['sequence', 'variant'] + angle_cols + 
                [col for col in combined_df.columns if col.startswith('Against')]]
feature_groups['mutations'] = mutation_cols
print(f"   Mutation features: {len(mutation_cols)}")

# Distance features (if available)
distance_cols = [col for col in combined_df.columns if col.startswith('residue_')]
if distance_cols:
    feature_groups['distances'] = distance_cols
    print(f"   Distance features: {len(distance_cols)}")
else:
    print("   No distance features found")

# Create feature matrix
all_features = []
feature_names = []

def clean_feature_name(name):
    """Clean feature name by removing special characters and spaces"""
    import re
    cleaned = re.sub(r'[^a-zA-Z0-9_]', '_', str(name))
    cleaned = re.sub(r'_+', '_', cleaned)
    cleaned = cleaned.strip('_')
    return cleaned

for group_name, cols in feature_groups.items():
    if cols:
        features = combined_df[cols].values.astype(float)
        all_features.append(features)
        # Create clean feature names
        for col in cols:
            clean_name = clean_feature_name(f"{group_name}_{col}")
            feature_names.append(clean_name)
        print(f"   Added {len(cols)} {group_name} features")

# Combine all features
X = np.hstack(all_features)
print(f"   Total features: {X.shape[1]}")

# ============================================================================
# 3. CREATE COMPREHENSIVE DENDROGRAM
# ============================================================================

print("\n4. Creating comprehensive dendrogram analysis...")

# Standardize the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Create figure with multiple subplots for different linkage methods
fig, axes = plt.subplots(2, 3, figsize=(20, 12))
fig.suptitle('Comprehensive Clustering Dendrogram Analysis: Angles + Mutations + Distances', 
             fontsize=16, y=0.98)

# Define linkage methods to test
linkage_methods = ['ward', 'complete', 'average', 'single', 'centroid']

for i, method in enumerate(linkage_methods):
    row = i // 3
    col = i % 3
    
    try:
        # Calculate distance matrix and linkage
        if method == 'ward':
            distance_matrix = pdist(X_scaled, metric='euclidean')
        else:
            distance_matrix = pdist(X_scaled, metric='euclidean')
        
        linkage_matrix = linkage(distance_matrix, method=method)
        
        # Create dendrogram
        dendrogram(linkage_matrix, 
                  labels=combined_df['sequence'].values,
                  ax=axes[row, col],
                  orientation='top',
                  leaf_rotation=90,
                  leaf_font_size=8)
        
        axes[row, col].set_title(f'Dendrogram - {method.capitalize()} Linkage')
        axes[row, col].set_xlabel('Samples')
        axes[row, col].set_ylabel('Distance')
        
    except Exception as e:
        print(f"   Error with {method} linkage: {e}")
        axes[row, col].text(0.5, 0.5, f'Error with {method}\nlinkage method', 
                           ha='center', va='center', transform=axes[row, col].transAxes)
        axes[row, col].set_title(f'Dendrogram - {method.capitalize()} Linkage (Error)')

plt.tight_layout()
plt.show()

# ============================================================================
# 4. CREATE DETAILED WARD DENDROGRAM
# ============================================================================

print("\n5. Creating detailed Ward dendrogram...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# Ward dendrogram with sample labels
distance_matrix = pdist(X_scaled, metric='euclidean')
linkage_matrix = linkage(distance_matrix, method='ward')

# Create dendrogram with sample labels
dendrogram(linkage_matrix, 
          labels=combined_df['sequence'].values,
          ax=ax1,
          orientation='top',
          leaf_rotation=90,
          leaf_font_size=8,
          color_threshold=0.7 * max(linkage_matrix[:, 2]))

ax1.set_title('Detailed Dendrogram - Ward Linkage (Sample Labels)', fontsize=14)
ax1.set_xlabel('Samples')
ax1.set_ylabel('Distance')

# Create dendrogram with variant labels
variant_labels = combined_df['variant'].values
dendrogram(linkage_matrix, 
          labels=variant_labels,
          ax=ax2,
          orientation='top',
          leaf_rotation=90,
          leaf_font_size=10,
          color_threshold=0.7 * max(linkage_matrix[:, 2]))

ax2.set_title('Detailed Dendrogram - Ward Linkage (Variant Labels)', fontsize=14)
ax2.set_xlabel('Variants')
ax2.set_ylabel('Distance')

plt.tight_layout()
plt.show()

# ============================================================================
# 5. ANALYZE CLUSTER STRUCTURE
# ============================================================================

print("\n6. Analyzing cluster structure...")

# Determine optimal number of clusters using elbow method
n_clusters_range = range(2, min(10, len(X_scaled) // 2))
inertias = []

for n_clusters in n_clusters_range:
    cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
    inertias.append(np.sum([np.sum((X_scaled[cluster_labels == i] - 
                                   np.mean(X_scaled[cluster_labels == i], axis=0))**2) 
                           for i in range(1, n_clusters + 1)]))

# Plot elbow curve
plt.figure(figsize=(10, 6))
plt.plot(n_clusters_range, inertias, 'bo-', linewidth=2, markersize=8)
plt.xlabel('Number of Clusters')
plt.ylabel('Inertia (Within-cluster Sum of Squares)')
plt.title('Elbow Method for Optimal Number of Clusters')
plt.grid(True, alpha=0.3)
plt.show()

# ============================================================================
# 6. CLUSTER DISTRIBUTION ANALYSIS
# ============================================================================

print("\n7. Cluster distribution analysis...")

fig, axes = plt.subplots(2, 2, figsize=(15, 12))

for i, n_clusters in enumerate([3, 4, 5, 6]):
    row = i // 2
    col = i % 2
    
    cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
    
    # Create cluster distribution plot
    cluster_variant_counts = pd.DataFrame({
        'variant': combined_df['variant'].values,
        'cluster': cluster_labels
    }).groupby(['cluster', 'variant']).size().unstack(fill_value=0)
    
    cluster_variant_counts.plot(kind='bar', ax=axes[row, col], colormap='tab10')
    axes[row, col].set_title(f'Cluster Distribution - {n_clusters} Clusters')
    axes[row, col].set_xlabel('Cluster')
    axes[row, col].set_ylabel('Number of Samples')
    axes[row, col].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    axes[row, col].tick_params(axis='x', rotation=0)

plt.tight_layout()
plt.show()

# ============================================================================
# 7. PRINT ANALYSIS SUMMARY
# ============================================================================

print("\n" + "=" * 70)
print("COMPREHENSIVE DENDROGRAM ANALYSIS SUMMARY")
print("=" * 70)

print(f"\nDataset Information:")
print(f"  Total samples: {len(combined_df)}")
print(f"  Total variants: {len(combined_df['variant'].unique())}")
print(f"  Variants: {', '.join(sorted(combined_df['variant'].unique()))}")

print(f"\nFeature Groups:")
for group_name, cols in feature_groups.items():
    if cols:
        print(f"  {group_name.capitalize()}: {len(cols)} features")

print(f"\nClustering Information:")
print(f"  Linkage method: Ward")
print(f"  Distance metric: Euclidean")
print(f"  Maximum distance: {max(linkage_matrix[:, 2]):.3f}")

# Calculate cophenetic correlation
distance_matrix = pdist(combined_df.drop(columns=['sequence', 'variant']).values)
cophenetic_distances = cophenet(linkage_matrix)
cophenetic_correlation = np.corrcoef(distance_matrix, cophenetic_distances)[0, 1]

print(f"  Cophenetic correlation: {cophenetic_correlation:.3f}")

print(f"\nSample Distribution by Variant:")
variant_counts = combined_df['variant'].value_counts()
for variant, count in variant_counts.items():
    print(f"  {variant}: {count} samples")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE!")
print("=" * 70) 