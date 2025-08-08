import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load and process the data
print("Loading mutation data...")
mutation = pd.read_csv('mutation/mutation-all.csv').drop(columns=['sequence']).groupby('variant').mean().reset_index()
mutation.to_csv('mutation/mutation-variant-count.csv', index=False)

# Set variant as index for heatmap
mutation_matrix = mutation.set_index('variant')

print(f"Data shape: {mutation_matrix.shape}")
print(f"Variants: {list(mutation_matrix.index)}")

# 1. Create full heatmap
plt.figure(figsize=(20, 8))
sns.heatmap(mutation_matrix, 
           annot=False,
           cmap='RdYlBu_r',  # Red for high frequency, blue for low
           vmin=0, vmax=1,
           cbar_kws={'label': 'Mutation Frequency'},
           xticklabels=True,
           yticklabels=True)

plt.title('COVID-19 Variant Mutation Frequency Heatmap', fontsize=16, pad=20)
plt.xlabel('Mutations', fontsize=12)
plt.ylabel('Variants', fontsize=12)
plt.xticks(rotation=90, fontsize=8)
plt.yticks(rotation=0, fontsize=10)
plt.tight_layout()
plt.savefig('variant_mutation_heatmap_full.png', dpi=300, bbox_inches='tight')
plt.close()

# 2. Create focused heatmap with high-variance mutations
print("Creating focused heatmap...")

# Calculate variance across variants for each mutation
mutation_variance = mutation_matrix.var(axis=0)
high_var_mutations = mutation_variance.nlargest(30)  # Top 30 most variable mutations

# Create focused heatmap
focused_data = mutation_matrix[high_var_mutations.index]

plt.figure(figsize=(15, 8))
sns.heatmap(focused_data, 
           annot=True,  # Show values since fewer columns
           fmt='.2f',
           cmap='RdYlBu_r',
           vmin=0, vmax=1,
           cbar_kws={'label': 'Mutation Frequency'},
           xticklabels=True,
           yticklabels=True)

plt.title('Top 30 Most Variable Mutations Across Variants', fontsize=16, pad=20)
plt.xlabel('Mutations', fontsize=12)
plt.ylabel('Variants', fontsize=12)
plt.xticks(rotation=90, fontsize=10)
plt.yticks(rotation=0, fontsize=12)
plt.tight_layout()
plt.savefig('variant_mutation_heatmap_focused.png', dpi=300, bbox_inches='tight')
plt.close()

# 3. Create variant similarity heatmap
print("Creating variant similarity heatmap...")

# Calculate correlation between variants
variant_corr = mutation_matrix.T.corr()

plt.figure(figsize=(8, 6))
sns.heatmap(variant_corr,
           annot=True,
           fmt='.3f',
           cmap='coolwarm',
           center=0,
           square=True,
           cbar_kws={'label': 'Correlation'})

plt.title('Variant Similarity (Mutation Pattern Correlation)', fontsize=14, pad=20)
plt.tight_layout()
plt.savefig('variant_similarity_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()

# Print summary statistics
print("\nSummary Statistics:")
print("=" * 50)
print(f"Total mutations analyzed: {mutation_matrix.shape[1]}")
print(f"Total variants: {mutation_matrix.shape[0]}")
print(f"Average mutation frequency: {mutation_matrix.values.mean():.3f}")
print(f"Most variable mutations:")
for i, (mut, var) in enumerate(high_var_mutations.head(10).items(), 1):
    print(f"  {i:2d}. {mut}: variance = {var:.3f}")

# 4. Create heatmap for residue region 319-541
print("Creating residue region 319-541 heatmap...")

# Filter mutations to only include those in the 319-541 region
# These would be mutations that affect residues in this range
residue_region_mutations = []
for col in mutation_matrix.columns:
    # Check if mutation name contains residue numbers in range 319-541
    # Common mutation format: Letter + Number + Letter (e.g., D614G, N856S)
    import re
    match = re.search(r'[A-Z](\d+)[A-Z-]', col)
    if match:
        residue_num = int(match.group(1))
        if 319 <= residue_num <= 541:
            residue_region_mutations.append(col)

print(f"Found {len(residue_region_mutations)} mutations in residue region 319-541")

if residue_region_mutations:
    # Create heatmap for this region
    region_data = mutation_matrix[residue_region_mutations]
    
    plt.figure(figsize=(16, 8))
    sns.heatmap(region_data, 
               annot=True,
               fmt='.2f',
               cmap='RdYlBu_r',
               vmin=0, vmax=1,
               cbar_kws={'label': 'Mutation Frequency'},
               xticklabels=True,
               yticklabels=True)

    plt.title('COVID-19 Variant Mutation Frequency Heatmap\n(Residue Region 319-541)', fontsize=16, pad=20)
    plt.xlabel('Mutations (Residue 319-541)', fontsize=12)
    plt.ylabel('Variants', fontsize=12)
    plt.xticks(rotation=90, fontsize=10)
    plt.yticks(rotation=0, fontsize=12)
    plt.tight_layout()
    plt.savefig('variant_mutation_heatmap_region_319_541.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print statistics for this region
    print(f"\nResidue Region 319-541 Statistics:")
    print(f"  • Mutations in region: {len(residue_region_mutations)}")
    print(f"  • Average frequency in region: {region_data.values.mean():.3f}")
    print(f"  • Most common mutations in region:")
    region_variance = region_data.var(axis=0)
    for i, (mut, var) in enumerate(region_variance.nlargest(10).items(), 1):
        avg_freq = region_data[mut].mean()
        print(f"    {i:2d}. {mut:<15} | variance={var:.3f} | avg_freq={avg_freq:.3f}")
else:
    print("No mutations found in residue region 319-541")

print("\nFiles created:")
print("- variant_mutation_heatmap_full.png (all mutations)")
print("- variant_mutation_heatmap_focused.png (top 30 variable mutations)")  
print("- variant_similarity_heatmap.png (variant correlations)")
if residue_region_mutations:
    print("- variant_mutation_heatmap_region_319_541.png (residue region 319-541)")
print("- mutation/mutation-variant-count.csv (processed data)")