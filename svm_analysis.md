# SVM Feature Interpretation Analysis

## Overview

**Total features available:** 240  
**Feature types:** 
- Mutation features: `['D614G', 'N856S', 'G1124V', 'T19R', 'T95I', 'G142D', 'E156G', 'F157-', 'R158-', 'R158G']`...
- Sliding window features: `sw_feature_*` (172 features derived from residue distances)

---

## Understanding SVM Coefficients and Importance Metrics

### What are SVM Coefficients?

In a **linear SVM**, the decision boundary is defined by the equation:
```
f(x) = w₁×x₁ + w₂×x₂ + ... + wₙ×xₙ + b = 0
```

Where:
- `wᵢ` = **coefficients** (weights) for each feature
- `xᵢ` = feature values  
- `b` = bias/intercept term

### Coefficient Interpretation

| Coefficient Sign | Meaning | Example from Our Data |
|-----------------|---------|----------------------|
| **Positive (+)** | Higher feature values push prediction toward **class +1** | `T20N` (+0.0438) favors positive class |
| **Negative (-)** | Higher feature values push prediction toward **class -1** | `N856S` (-0.4125) strongly favors negative class |
| **Large magnitude** | Stronger influence on classification | `N856S` (-0.4125) >> `T20N` (+0.0438) |

### Importance Metric

**Formula:** `importance = abs(coefficient)`

The **importance** represents:
- How much that feature influences the decision boundary
- Independent of which class it favors  
- Allows ranking features by their overall predictive power

### Real Example: Omicron-BA.1 Binding Affinity Classification

| Feature | Coefficient | Importance | Biological Interpretation |
|---------|-------------|------------|-------------------------|
| N856S | -0.4125 | 0.4125 | **Very strong** predictor; presence indicates LOWER binding affinity than Omicron-BA.1 |
| D614G | -0.3876 | 0.3876 | **Strong** predictor; presence indicates LOWER binding affinity than Omicron-BA.1 |
| V1176F | +0.0209 | 0.0209 | **Weak** predictor; presence slightly indicates HIGHER binding affinity than Omicron-BA.1 |
| sw_feature_399_min | +0.0195 | 0.0195 | Structural pattern; higher minimum distances indicate HIGHER binding affinity than Omicron-BA.1 |

### Decision Process

For any new sequence, the SVM calculates:
```
Score = (-0.4125 × N856S_value) + (-0.3876 × D614G_value) + 
        (+0.0209 × V1176F_value) + ... + bias
```

**Classification rule:**
- **Score > 0** → **HIGHER** binding affinity than the reference variant
- **Score < 0** → **LOWER** binding affinity than the reference variant

### Key Insights from Our Analysis

1. **Mutation features dominate:** Top predictors are sequence mutations (N856S, D614G, R346K)
2. **Sliding window features contribute:** Structural distance patterns add predictive value
3. **Feature consistency:** Features appearing in 4/4 models are robust cross-variant biomarkers
4. **Biological relevance:** High-importance features correspond to known variant-defining mutations
5. **Binding affinity prediction:** Models predict whether sequences have higher or lower binding affinity than reference variants

---

## Feature Importance Analysis by Binding Affinity Classification

### Model 1: Wuhan-A Binding Affinity Classification

#### Top 15 Most Important Features

| Rank | Feature | Coefficient | Importance |
|------|---------|-------------|------------|
| 1 | T20N | 0.0438 | 0.0438 |
| 2 | V1176F | 0.0438 | 0.0438 |
| 3 | T1027I | 0.0416 | 0.0416 |
| 4 | R190S | 0.0416 | 0.0416 |
| 5 | D138Y | 0.0416 | 0.0416 |
| 6 | P26S | 0.0416 | 0.0416 |
| 7 | L18F | 0.0416 | 0.0416 |
| 8 | K417T | 0.0398 | 0.0398 |
| 9 | H655Y | 0.0299 | 0.0299 |
| 10 | sw_feature_499_std | 0.0283 | 0.0283 |
| 11 | sw_feature_504_std | 0.0282 | 0.0282 |
| 12 | sw_feature_504_max | 0.0281 | 0.0281 |
| 13 | sw_feature_519_max | -0.0280 | 0.0280 |
| 14 | N856S | -0.0269 | 0.0269 |
| 15 | sw_feature_434_std | 0.0253 | 0.0253 |

#### Interpretation:
- **Strongest positive influence:** T20N (indicates HIGHER binding affinity than Wuhan-A)
- **Features favoring LOWER binding affinity than Wuhan-A (-1):** 2 features
- **Strongest negative influence:** sw_feature_519_max (indicates LOWER binding affinity than Wuhan-A)

---

### Model 2: Beta-B.1.351 Binding Affinity Classification

#### Top 15 Most Important Features

| Rank | Feature | Coefficient | Importance |
|------|---------|-------------|------------|
| 1 | R346K | -0.1068 | 0.1068 |
| 2 | sw_feature_449_std | -0.0759 | 0.0759 |
| 3 | sw_feature_449_max | -0.0731 | 0.0731 |
| 4 | sw_feature_329_std | 0.0659 | 0.0659 |
| 5 | sw_feature_469_std | 0.0657 | 0.0657 |
| 6 | sw_feature_329_max | 0.0622 | 0.0622 |
| 7 | sw_feature_324_max | 0.0621 | 0.0621 |
| 8 | sw_feature_364_min | 0.0545 | 0.0545 |
| 9 | sw_feature_329_mean | 0.0544 | 0.0544 |
| 10 | sw_feature_324_min | 0.0533 | 0.0533 |
| 11 | K417N | -0.0491 | 0.0491 |
| 12 | sw_feature_469_max | 0.0466 | 0.0466 |
| 13 | sw_feature_474_mean | 0.0464 | 0.0464 |
| 14 | G446S | -0.0462 | 0.0462 |
| 15 | N440K | -0.0462 | 0.0462 |

#### Interpretation:
- **Strongest positive influence:** sw_feature_329_std (indicates HIGHER binding affinity than Beta-B.1.351)
- **Features favoring LOWER binding affinity than Beta-B.1.351 (-1):** 6 features
- **Strongest negative influence:** R346K (indicates LOWER binding affinity than Beta-B.1.351)

---

### Model 3: Omicron-BA.1 Binding Affinity Classification

#### Top 15 Most Important Features

| Rank | Feature | Coefficient | Importance |
|------|---------|-------------|------------|
| 1 | N856S | -0.4125 | 0.4125 |
| 2 | D614G | -0.3876 | 0.3876 |
| 3 | V1176F | 0.0209 | 0.0209 |
| 4 | T20N | 0.0209 | 0.0209 |
| 5 | L452R | 0.0203 | 0.0203 |
| 6 | P681R | 0.0203 | 0.0203 |
| 7 | D950N | 0.0203 | 0.0203 |
| 8 | D138Y | 0.0199 | 0.0199 |
| 9 | T1027I | 0.0199 | 0.0199 |
| 10 | L18F | 0.0199 | 0.0199 |
| 11 | P26S | 0.0199 | 0.0199 |
| 12 | R190S | 0.0199 | 0.0199 |
| 13 | sw_feature_399_min | 0.0195 | 0.0195 |
| 14 | sw_feature_394_min | 0.0195 | 0.0195 |
| 15 | sw_feature_384_std | 0.0193 | 0.0193 |

#### Interpretation:
- **Strongest positive influence:** V1176F (indicates HIGHER binding affinity than Omicron-BA.1)
- **Features favoring LOWER binding affinity than Omicron-BA.1 (-1):** 2 features
- **Strongest negative influence:** N856S (indicates LOWER binding affinity than Omicron-BA.1)

---

### Model 4: Delta-B.1.617.2 Binding Affinity Classification

#### Top 15 Most Important Features

| Rank | Feature | Coefficient | Importance |
|------|---------|-------------|------------|
| 1 | N856S | -0.4049 | 0.4049 |
| 2 | D614G | -0.3949 | 0.3949 |
| 3 | T20N | 0.0234 | 0.0234 |
| 4 | V1176F | 0.0234 | 0.0234 |
| 5 | T1027I | 0.0222 | 0.0222 |
| 6 | L18F | 0.0222 | 0.0222 |
| 7 | P26S | 0.0222 | 0.0222 |
| 8 | D138Y | 0.0222 | 0.0222 |
| 9 | R190S | 0.0222 | 0.0222 |
| 10 | K417T | 0.0213 | 0.0213 |
| 11 | sw_feature_489_min | 0.0166 | 0.0166 |
| 12 | sw_feature_484_max | 0.0154 | 0.0154 |
| 13 | sw_feature_494_min | 0.0149 | 0.0149 |
| 14 | sw_feature_439_mean | 0.0140 | 0.0140 |
| 15 | sw_feature_484_min | 0.0139 | 0.0139 |

#### Interpretation:
- **Strongest positive influence:** T20N (indicates HIGHER binding affinity than Delta-B.1.617.2)
- **Features favoring LOWER binding affinity than Delta-B.1.617.2 (-1):** 2 features
- **Strongest negative influence:** N856S (indicates LOWER binding affinity than Delta-B.1.617.2)

---

## Summary: Cross-Model Feature Importance

### Top 20 Features (Average Importance Across 4 Binding Affinity Classification Models)

| Rank | Feature | Avg Importance | Frequency |
|------|---------|----------------|-----------|
| 1 | N856S | 0.2111 | 4/4 |
| 2 | D614G | 0.1976 | 4/4 |
| 3 | R346K | 0.0313 | 4/4 |
| 4 | sw_feature_469_std | 0.0240 | 4/4 |
| 5 | sw_feature_449_std | 0.0224 | 4/4 |
| 6 | T20N | 0.0220 | 4/4 |
| 7 | V1176F | 0.0220 | 4/4 |
| 8 | sw_feature_449_max | 0.0215 | 4/4 |
| 9 | sw_feature_324_max | 0.0212 | 4/4 |
| 10 | T1027I | 0.0209 | 4/4 |
| 11 | R190S | 0.0209 | 4/4 |
| 12 | D138Y | 0.0209 | 4/4 |
| 13 | P26S | 0.0209 | 4/4 |
| 14 | L18F | 0.0209 | 4/4 |
| 15 | sw_feature_469_max | 0.0205 | 4/4 |
| 16 | K417T | 0.0196 | 4/4 |
| 17 | sw_feature_474_mean | 0.0191 | 4/4 |
| 18 | sw_feature_324_min | 0.0185 | 4/4 |
| 19 | sw_feature_329_std | 0.0178 | 4/4 |
| 20 | sw_feature_364_min | 0.0174 | 4/4 |

### Key Findings

**Most consistently important features across all binding affinity classification models:**

1. **N856S** - avg importance = 0.2111 (appears in 4/4 models) - Strong predictor of binding affinity differences
2. **D614G** - avg importance = 0.1976 (appears in 4/4 models) - Strong predictor of binding affinity differences
3. **R346K** - avg importance = 0.0313 (appears in 4/4 models) - Moderate predictor of binding affinity differences
4. **sw_feature_469_std** - avg importance = 0.0240 (appears in 4/4 models) - Structural predictor of binding affinity differences
5. **sw_feature_449_std** - avg importance = 0.0224 (appears in 4/4 models) - Structural predictor of binding affinity differences

### Insights

- **Mutation features dominate:** The top 3 features are all sequence mutations (N856S, D614G, R346K) that strongly predict binding affinity differences
- **Sliding window features contribute:** Multiple sliding window features appear in top rankings, indicating that structural distance patterns predict binding affinity differences
- **Consistent patterns:** Features like T20N, V1176F, and several sliding window statistics appear consistently across binding affinity classification models
- **Structural regions of interest:** Features in the 320-470 range (sw_feature_324_*, sw_feature_329_*, sw_feature_449_*, sw_feature_469_*) show high importance for binding affinity prediction
- **Binding affinity classification:** Models successfully distinguish between sequences with higher vs. lower binding affinity than reference variants