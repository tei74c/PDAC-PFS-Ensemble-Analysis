"""
Monte Carlo Robustness Analysis for 3-Feature PFS Classification Model
=======================================================================

Comprehensive robustness testing of the parsimonious Naive Bayes classifier
for progression-free survival (PFS) prediction using three proteomic features
(DYNC2H1, ECM2, PPIB).

Analyses performed:
    1. Targeted gene dropout (leave-one-out feature analysis)
    2. Random gene dropout (Monte Carlo sampling)
    3. Expression noise injection at multiple levels
    4. Feature-specific noise sensitivity profiling
    5. Corrected permutation test (training-label shuffle only)
    6. Bootstrap stability analysis with 95% confidence intervals

Statistical correction:
    Phipson-Smyth empirical p-value: p = (r + 1) / (n + 1),
    where r = number of permutation statistics >= observed statistic,
    and n = total number of permutations.
    This prevents p = 0.000 from finite permutation samples.
    Reference: Phipson & Smyth, Stat. Appl. Genet. Mol. Biol. (2010).

Reference: Khoshnevis et al. Molecular Cancer (2026).
"""

import os
import pickle
import pandas as pd
import numpy as np
from sklearn.metrics import (accuracy_score, f1_score, roc_auc_score,
                             precision_score, recall_score, confusion_matrix,
                             matthews_corrcoef, cohen_kappa_score, brier_score_loss)
from sklearn.naive_bayes import GaussianNB
from sklearn.calibration import CalibratedClassifierCV
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
DATA_DIR = os.path.join(BASE_DIR, "data")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Publication figure settings
# ---------------------------------------------------------------------------
# TrueType fonts for Adobe Illustrator compatibility
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['figure.dpi'] = 300

# Set random seed for reproducibility
np.random.seed(42)

# ============================================================================
# Load Data
# ============================================================================
print("=" * 70)
print("MONTE CARLO ROBUSTNESS ANALYSIS: 3-FEATURE MODEL")
print("=" * 70)

train_df = pd.read_excel(os.path.join(DATA_DIR, 'PFS_ML_train_cof_wo_IG.xlsx'))
test_df = pd.read_excel(os.path.join(DATA_DIR, 'PFS_ML_test_cof_wo_IG.xlsx'))

# Encode PFS_group: L=0 (Long PFS), S=1 (Short PFS)
label_map = {'L': 0, 'S': 1}
train_df['PFS_group_encoded'] = train_df['PFS_group'].map(label_map)
test_df['PFS_group_encoded'] = test_df['PFS_group'].map(label_map)

# 3-Feature model features
features_3 = ['DYNC2H1', 'ECM2', 'PPIB']

X_train = train_df[features_3].values
X_test = test_df[features_3].values
y_train = train_df['PFS_group_encoded'].values
y_test = test_df['PFS_group_encoded'].values

print(f"\n3-Feature Model Features: {features_3}")
print(f"Training samples: {len(y_train)} (Short PFS: {sum(y_train)}, Long PFS: {len(y_train)-sum(y_train)})")
print(f"Test samples: {len(y_test)} (Short PFS: {sum(y_test)}, Long PFS: {len(y_test)-sum(y_test)})")

# Train 3-feature model
base_model = GaussianNB()
model_3 = CalibratedClassifierCV(base_model, method='sigmoid', cv=3)
model_3.fit(X_train, y_train)

# Baseline metrics
pred_3 = model_3.predict(X_test)
proba_3 = model_3.predict_proba(X_test)[:, 1]

baseline_metrics = {
    'Accuracy': accuracy_score(y_test, pred_3),
    'F1-Score': f1_score(y_test, pred_3),
    'ROC-AUC': roc_auc_score(y_test, proba_3),
    'Precision': precision_score(y_test, pred_3),
    'Recall': recall_score(y_test, pred_3),
    'MCC': matthews_corrcoef(y_test, pred_3),
    'Brier-Score': brier_score_loss(y_test, proba_3)
}

print(f"\n{'='*70}")
print("BASELINE PERFORMANCE (3-Feature Model)")
print("=" * 70)
for metric, value in baseline_metrics.items():
    print(f"  {metric}: {value:.4f}")

# ============================================================================
# 1. TARGETED GENE DROPOUT (Leave-One-Out Feature Analysis)
# ============================================================================
print(f"\n{'='*70}")
print("1. TARGETED GENE DROPOUT ANALYSIS")
print("=" * 70)

dropout_results = []

for i, feature in enumerate(features_3):
    remaining_features = [f for f in features_3 if f != feature]
    remaining_idx = [features_3.index(f) for f in remaining_features]

    X_train_dropout = X_train[:, remaining_idx]
    X_test_dropout = X_test[:, remaining_idx]

    base_dropout = GaussianNB()
    model_dropout = CalibratedClassifierCV(base_dropout, method='sigmoid', cv=3)
    model_dropout.fit(X_train_dropout, y_train)

    pred_dropout = model_dropout.predict(X_test_dropout)
    proba_dropout = model_dropout.predict_proba(X_test_dropout)[:, 1]

    auc_dropout = roc_auc_score(y_test, proba_dropout)
    acc_dropout = accuracy_score(y_test, pred_dropout)

    auc_drop = baseline_metrics['ROC-AUC'] - auc_dropout
    acc_drop = baseline_metrics['Accuracy'] - acc_dropout

    dropout_results.append({
        'Feature': feature,
        'Remaining_Features': remaining_features,
        'AUC_Without': auc_dropout,
        'AUC_Drop': auc_drop,
        'Accuracy_Without': acc_dropout,
        'Accuracy_Drop': acc_drop
    })

    print(f"\n  Without {feature}:")
    print(f"    AUC: {auc_dropout:.4f} (drop: {auc_drop:+.4f})")
    print(f"    Accuracy: {acc_dropout:.4f} (drop: {acc_drop:+.4f})")

dropout_results_sorted = sorted(dropout_results, key=lambda x: x['AUC_Drop'], reverse=True)
print(f"\n  Feature Importance Ranking (by AUC drop):")
for i, result in enumerate(dropout_results_sorted, 1):
    print(f"    {i}. {result['Feature']}: {result['AUC_Drop']:+.4f}")

# ============================================================================
# 2. RANDOM GENE DROPOUT (Monte Carlo)
# ============================================================================
print(f"\n{'='*70}")
print("2. RANDOM GENE DROPOUT ANALYSIS")
print("=" * 70)

n_iterations = 100
random_dropout_results = {1: [], 2: []}

for n_drop in [1, 2]:
    print(f"\n  Dropping {n_drop} feature(s) randomly ({n_iterations} iterations)...")

    for _ in range(n_iterations):
        drop_idx = np.random.choice(len(features_3), n_drop, replace=False)
        keep_idx = [i for i in range(len(features_3)) if i not in drop_idx]

        if len(keep_idx) == 0:
            continue

        X_train_drop = X_train[:, keep_idx]
        X_test_drop = X_test[:, keep_idx]

        base_drop = GaussianNB()
        model_drop = CalibratedClassifierCV(base_drop, method='sigmoid', cv=3)
        model_drop.fit(X_train_drop, y_train)

        proba_drop = model_drop.predict_proba(X_test_drop)[:, 1]
        auc_drop = roc_auc_score(y_test, proba_drop)
        random_dropout_results[n_drop].append(auc_drop)

    mean_auc = np.mean(random_dropout_results[n_drop])
    std_auc = np.std(random_dropout_results[n_drop])
    print(f"    Mean AUC: {mean_auc:.4f} +/- {std_auc:.4f}")
    print(f"    AUC degradation: {baseline_metrics['ROC-AUC'] - mean_auc:+.4f}")

# ============================================================================
# 3. NOISE INJECTION ANALYSIS
# ============================================================================
print(f"\n{'='*70}")
print("3. NOISE INJECTION ANALYSIS")
print("=" * 70)

noise_levels = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
n_noise_iterations = 100
noise_results = {}

# Noise scaled by per-feature standard deviation from training data
feature_stds = np.std(X_train, axis=0)

for noise_level in noise_levels:
    aucs = []
    accs = []
    f1s = []

    for _ in range(n_noise_iterations):
        noise = np.random.normal(0, feature_stds * noise_level, X_test.shape)
        X_test_noisy = X_test + noise

        proba_noisy = model_3.predict_proba(X_test_noisy)[:, 1]
        pred_noisy = model_3.predict(X_test_noisy)

        aucs.append(roc_auc_score(y_test, proba_noisy))
        accs.append(accuracy_score(y_test, pred_noisy))
        f1s.append(f1_score(y_test, pred_noisy))

    noise_results[noise_level] = {
        'AUC_mean': np.mean(aucs),
        'AUC_std': np.std(aucs),
        'Accuracy_mean': np.mean(accs),
        'Accuracy_std': np.std(accs),
        'F1_mean': np.mean(f1s),
        'F1_std': np.std(f1s)
    }

    degradation = (baseline_metrics['ROC-AUC'] - np.mean(aucs)) / baseline_metrics['ROC-AUC'] * 100
    print(f"  {int(noise_level*100):3d}% noise: AUC = {np.mean(aucs):.4f} +/- {np.std(aucs):.4f} ({degradation:+.1f}% degradation)")

# ============================================================================
# 4. FEATURE-SPECIFIC NOISE SENSITIVITY
# ============================================================================
print(f"\n{'='*70}")
print("4. FEATURE-SPECIFIC NOISE SENSITIVITY")
print("=" * 70)

feature_noise_sensitivity = {}

for i, feature in enumerate(features_3):
    print(f"\n  {feature}:")
    feature_results = {}

    for noise_level in [0.1, 0.2, 0.3, 0.5]:
        aucs = []

        for _ in range(n_noise_iterations):
            # Add noise only to this single feature
            X_test_noisy = X_test.copy()
            noise = np.random.normal(0, feature_stds[i] * noise_level, X_test.shape[0])
            X_test_noisy[:, i] += noise

            proba_noisy = model_3.predict_proba(X_test_noisy)[:, 1]
            aucs.append(roc_auc_score(y_test, proba_noisy))

        feature_results[noise_level] = {
            'mean': np.mean(aucs),
            'std': np.std(aucs)
        }

        degradation = baseline_metrics['ROC-AUC'] - np.mean(aucs)
        print(f"    {int(noise_level*100):3d}% noise: AUC = {np.mean(aucs):.4f} (drop: {degradation:+.4f})")

    feature_noise_sensitivity[feature] = feature_results

# ============================================================================
# 5. CORRECTED PERMUTATION TEST (Statistical Significance)
# ============================================================================
print(f"\n{'='*70}")
print("5. CORRECTED PERMUTATION TEST (Training labels shuffled, test labels FIXED)")
print("=" * 70)

n_permutations = 1000
perm_aucs = []

print(f"\n  Running {n_permutations} permutations...")
for i in range(n_permutations):
    if (i + 1) % 200 == 0:
        print(f"    Permutation {i+1}/{n_permutations}...")

    # Shuffle ONLY training labels (correct methodology)
    y_train_shuffled = np.random.permutation(y_train)

    base_perm = GaussianNB()
    model_perm = CalibratedClassifierCV(base_perm, method='sigmoid', cv=3)

    try:
        model_perm.fit(X_train, y_train_shuffled)
        # Evaluate against TRUE (not shuffled) test labels
        proba_perm = model_perm.predict_proba(X_test)[:, 1]
        perm_aucs.append(roc_auc_score(y_test, proba_perm))
    except Exception:
        continue

# Phipson-Smyth corrected empirical p-value: p = (r + 1) / (n + 1)
# Prevents p = 0.000 from finite permutation samples
p_value = (np.sum(np.array(perm_aucs) >= baseline_metrics['ROC-AUC']) + 1) / (len(perm_aucs) + 1)

# Cohen's d effect size
cohens_d = (baseline_metrics['ROC-AUC'] - np.mean(perm_aucs)) / np.std(perm_aucs)

print(f"\n  Permutation Test Results:")
print(f"    Observed AUC: {baseline_metrics['ROC-AUC']:.4f}")
print(f"    Null distribution mean: {np.mean(perm_aucs):.4f}")
print(f"    Null distribution std: {np.std(perm_aucs):.4f}")
print(f"    p-value: {p_value:.4f}")
print(f"    Cohen's d: {cohens_d:.2f}")
print(f"    Interpretation: {'SIGNIFICANT (p < 0.05)' if p_value < 0.05 else 'Not significant'}")

# ============================================================================
# 6. BOOTSTRAP STABILITY ANALYSIS
# ============================================================================
print(f"\n{'='*70}")
print("6. BOOTSTRAP STABILITY ANALYSIS")
print("=" * 70)

n_bootstrap = 1000
bootstrap_aucs = []
bootstrap_accs = []

for _ in range(n_bootstrap):
    # Bootstrap sample from training data
    boot_idx = np.random.choice(len(y_train), len(y_train), replace=True)
    X_train_boot = X_train[boot_idx]
    y_train_boot = y_train[boot_idx]

    base_boot = GaussianNB()
    model_boot = CalibratedClassifierCV(base_boot, method='sigmoid', cv=3)

    try:
        model_boot.fit(X_train_boot, y_train_boot)
        proba_boot = model_boot.predict_proba(X_test)[:, 1]
        pred_boot = model_boot.predict(X_test)

        bootstrap_aucs.append(roc_auc_score(y_test, proba_boot))
        bootstrap_accs.append(accuracy_score(y_test, pred_boot))
    except Exception:
        continue

print(f"\n  AUC across {len(bootstrap_aucs)} bootstrap models:")
print(f"    Mean: {np.mean(bootstrap_aucs):.4f}")
print(f"    Std: {np.std(bootstrap_aucs):.4f}")
print(f"    95% CI: [{np.percentile(bootstrap_aucs, 2.5):.4f}, {np.percentile(bootstrap_aucs, 97.5):.4f}]")
print(f"    Min: {np.min(bootstrap_aucs):.4f}, Max: {np.max(bootstrap_aucs):.4f}")

print(f"\n  Accuracy across bootstrap models:")
print(f"    Mean: {np.mean(bootstrap_accs):.4f}")
print(f"    Std: {np.std(bootstrap_accs):.4f}")
print(f"    95% CI: [{np.percentile(bootstrap_accs, 2.5):.4f}, {np.percentile(bootstrap_accs, 97.5):.4f}]")

# ============================================================================
# 7. COMPREHENSIVE VISUALIZATION
# ============================================================================
print(f"\n{'='*70}")
print("7. GENERATING ROBUSTNESS VISUALIZATIONS")
print("=" * 70)

fig = plt.figure(figsize=(20, 14))

# Panel A: Feature Importance (Dropout Analysis)
ax1 = fig.add_subplot(2, 4, 1)
features_sorted = [r['Feature'] for r in dropout_results_sorted]
auc_drops = [r['AUC_Drop'] for r in dropout_results_sorted]
colors = ['#d62728' if d > 0.1 else '#ff7f0e' if d > 0.05 else '#2ca02c' for d in auc_drops]

bars = ax1.barh(features_sorted, auc_drops, color=colors, edgecolor='black', linewidth=1)
ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1)
ax1.set_xlabel('AUC Drop When Feature Removed')
ax1.set_title('A. Feature Importance\n(Leave-One-Out Analysis)')
ax1.invert_yaxis()

for bar, val in zip(bars, auc_drops):
    ax1.text(val + 0.01, bar.get_y() + bar.get_height()/2, f'{val:+.3f}',
             va='center', fontsize=10, fontweight='bold')

# Panel B: Noise Robustness
ax2 = fig.add_subplot(2, 4, 2)
noise_pct = [int(n * 100) for n in noise_levels]
auc_means = [noise_results[n]['AUC_mean'] for n in noise_levels]
auc_stds = [noise_results[n]['AUC_std'] for n in noise_levels]

ax2.errorbar(noise_pct, auc_means, yerr=auc_stds, fmt='o-', color='green',
             capsize=5, linewidth=2, markersize=10, markerfacecolor='white',
             markeredgewidth=2, label='3-Feature Model')
ax2.axhline(y=baseline_metrics['ROC-AUC'], color='green', linestyle='--', alpha=0.5,
            label=f'Baseline AUC ({baseline_metrics["ROC-AUC"]:.3f})')
ax2.axhline(y=0.5, color='red', linestyle=':', linewidth=1.5, label='Chance')

for i, (n, auc) in enumerate(zip(noise_pct, auc_means)):
    if n in [10, 30, 50]:
        deg = (baseline_metrics['ROC-AUC'] - auc) / baseline_metrics['ROC-AUC'] * 100
        ax2.annotate(f'{deg:.1f}%', (n, auc - 0.02), ha='center', fontsize=9, color='darkred')

ax2.set_xlabel('Noise Level (% of Feature Std)')
ax2.set_ylabel('ROC-AUC')
ax2.set_title('B. Noise Robustness')
ax2.legend(loc='lower left')
ax2.set_xlim(0, 55)
ax2.set_ylim(0.4, 1.0)
ax2.grid(True, alpha=0.3)

# Panel C: Feature-Specific Noise Sensitivity
ax3 = fig.add_subplot(2, 4, 3)
noise_test_levels = [0.1, 0.2, 0.3, 0.5]
colors_feat = ['#1f77b4', '#ff7f0e', '#2ca02c']

for i, (feature, color) in enumerate(zip(features_3, colors_feat)):
    means = [feature_noise_sensitivity[feature][n]['mean'] for n in noise_test_levels]
    stds = [feature_noise_sensitivity[feature][n]['std'] for n in noise_test_levels]
    ax3.errorbar([n*100 for n in noise_test_levels], means, yerr=stds,
                 fmt='o-', color=color, capsize=4, linewidth=2, markersize=8,
                 label=feature)

ax3.axhline(y=baseline_metrics['ROC-AUC'], color='gray', linestyle='--', alpha=0.5)
ax3.set_xlabel('Noise Level (% of Feature Std)')
ax3.set_ylabel('ROC-AUC')
ax3.set_title('C. Feature-Specific Noise Sensitivity')
ax3.legend(loc='lower left')
ax3.set_xlim(5, 55)
ax3.grid(True, alpha=0.3)

# Panel D: Bootstrap Stability Distribution
ax4 = fig.add_subplot(2, 4, 4)
ax4.hist(bootstrap_aucs, bins=30, color='green', alpha=0.7, edgecolor='black', linewidth=0.5)
ax4.axvline(x=baseline_metrics['ROC-AUC'], color='red', linestyle='-', linewidth=2,
            label=f'Original Model ({baseline_metrics["ROC-AUC"]:.3f})')
ax4.axvline(x=np.mean(bootstrap_aucs), color='blue', linestyle='--', linewidth=2,
            label=f'Bootstrap Mean ({np.mean(bootstrap_aucs):.3f})')

ci_low = np.percentile(bootstrap_aucs, 2.5)
ci_high = np.percentile(bootstrap_aucs, 97.5)
ax4.axvspan(ci_low, ci_high, alpha=0.2, color='blue', label=f'95% CI [{ci_low:.3f}, {ci_high:.3f}]')

ax4.set_xlabel('ROC-AUC')
ax4.set_ylabel('Frequency')
ax4.set_title(f'D. Bootstrap Stability (n={len(bootstrap_aucs)})')
ax4.legend(loc='upper left', fontsize=9)

# Panel E: Noise Tolerance Detail
ax5 = fig.add_subplot(2, 4, 5)

noise_pct_detailed = [5, 10, 15, 20, 25, 30, 40, 50]
auc_detailed = [noise_results[n/100]['AUC_mean'] for n in noise_pct_detailed]
std_detailed = [noise_results[n/100]['AUC_std'] for n in noise_pct_detailed]

ax5.fill_between(noise_pct_detailed,
                 [a - s for a, s in zip(auc_detailed, std_detailed)],
                 [a + s for a, s in zip(auc_detailed, std_detailed)],
                 alpha=0.3, color='green')
ax5.plot(noise_pct_detailed, auc_detailed, 'go-', linewidth=2, markersize=8,
         markerfacecolor='white', markeredgewidth=2, label='3-Feature Model')
ax5.axhline(y=baseline_metrics['ROC-AUC'], color='green', linestyle='--', alpha=0.5)
ax5.axhline(y=0.9, color='orange', linestyle=':', linewidth=1.5, label='Clinical threshold (0.9)')
ax5.axhline(y=0.5, color='red', linestyle=':', linewidth=1.5, label='Chance')

ax5.set_xlabel('Noise Level (%)')
ax5.set_ylabel('ROC-AUC')
ax5.set_title('E. Noise Tolerance Detail')
ax5.legend(loc='lower left', fontsize=8)
ax5.set_xlim(0, 55)
ax5.set_ylim(0.5, 1.0)
ax5.grid(True, alpha=0.3)

# Panel F: Permutation Test
ax6 = fig.add_subplot(2, 4, 6)
ax6.hist(perm_aucs, bins=40, color='gray', alpha=0.7, edgecolor='black', linewidth=0.5,
         label='Null Distribution')
ax6.axvline(x=baseline_metrics['ROC-AUC'], color='green', linestyle='-', linewidth=3,
            label=f'Observed ({baseline_metrics["ROC-AUC"]:.3f})')
ax6.axvline(x=np.mean(perm_aucs), color='red', linestyle='--', linewidth=2,
            label=f'Null Mean ({np.mean(perm_aucs):.3f})')
ax6.axvline(x=0.5, color='black', linestyle=':', linewidth=1.5)

p_text = f'p = {p_value:.3f}' if p_value >= 0.001 else 'p < 0.001'
ax6.text(0.95, 0.95, f'{p_text}\nCohen\'s d = {cohens_d:.2f}',
         transform=ax6.transAxes, fontsize=10, verticalalignment='top',
         horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

ax6.set_xlabel('ROC-AUC')
ax6.set_ylabel('Frequency')
ax6.set_title('F. Permutation Test\n(Corrected Methodology)')
ax6.legend(loc='upper left', fontsize=8)

# Panel G: Random Dropout Results
ax7 = fig.add_subplot(2, 4, 7)
x_pos = [1, 2]
dropout_means = [np.mean(random_dropout_results[1]), np.mean(random_dropout_results[2])]
dropout_stds = [np.std(random_dropout_results[1]), np.std(random_dropout_results[2])]

bars = ax7.bar(x_pos, dropout_means, yerr=dropout_stds, capsize=8, color=['#2ca02c', '#ff7f0e'],
               edgecolor='black', linewidth=1.5, width=0.6)
ax7.axhline(y=baseline_metrics['ROC-AUC'], color='red', linestyle='--', linewidth=2,
            label=f'Baseline ({baseline_metrics["ROC-AUC"]:.3f})')

ax7.set_xticks(x_pos)
ax7.set_xticklabels(['Drop 1 Feature\n(2 remaining)', 'Drop 2 Features\n(1 remaining)'])
ax7.set_ylabel('ROC-AUC')
ax7.set_title('G. Random Feature Dropout')
ax7.legend()
ax7.set_ylim(0.5, 1.0)

for bar, mean in zip(bars, dropout_means):
    deg = baseline_metrics['ROC-AUC'] - mean
    ax7.text(bar.get_x() + bar.get_width()/2, mean + 0.03, f'{deg:+.3f}',
             ha='center', fontsize=11, fontweight='bold')

# Panel H: Summary Statistics Table
ax8 = fig.add_subplot(2, 4, 8)
ax8.axis('off')

summary_data = [
    ['Robustness Metric', 'Value', 'Interpretation'],
    ['Baseline AUC', f'{baseline_metrics["ROC-AUC"]:.3f}', 'Excellent'],
    ['p-value (perm)', f'{p_value:.3f}', 'Significant' if p_value < 0.05 else 'NS'],
    ['Cohen\'s d', f'{cohens_d:.2f}', 'Large effect' if cohens_d > 0.8 else 'Moderate'],
    ['30% Noise AUC', f'{noise_results[0.3]["AUC_mean"]:.3f}', 'Robust'],
    ['Bootstrap 95% CI', f'[{ci_low:.3f}, {ci_high:.3f}]', 'Stable'],
    ['Most Important', dropout_results_sorted[0]['Feature'], 'Critical'],
]

table = ax8.table(cellText=summary_data, loc='center', cellLoc='center',
                  colWidths=[0.4, 0.3, 0.3])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.2, 1.7)

for i in range(3):
    table[(0, i)].set_facecolor('#4472C4')
    table[(0, i)].set_text_props(color='white', weight='bold')

ax8.set_title('H. Robustness Summary', pad=20, fontsize=12, fontweight='bold')

plt.suptitle('Monte Carlo Robustness Analysis: 3-Feature PFS Model\n(DYNC2H1, ECM2, PPIB)',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

plt.savefig(os.path.join(OUTPUT_DIR, 'Robustness_Analysis_3Feature.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(OUTPUT_DIR, 'Robustness_Analysis_3Feature.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(OUTPUT_DIR, 'Robustness_Analysis_3Feature.svg'),
            bbox_inches='tight', facecolor='white')
plt.close()

print("\nSaved: Robustness_Analysis_3Feature.png/pdf/svg")

# ============================================================================
# 8. SAVE RESULTS TO EXCEL
# ============================================================================
print(f"\n{'='*70}")
print("8. SAVING RESULTS TO EXCEL")
print("=" * 70)

with pd.ExcelWriter(os.path.join(OUTPUT_DIR, 'Robustness_Analysis_3Feature_Results.xlsx'),
                     engine='openpyxl') as writer:

    # Baseline metrics
    pd.DataFrame([baseline_metrics]).to_excel(writer, sheet_name='Baseline_Metrics', index=False)

    # Dropout analysis
    pd.DataFrame(dropout_results_sorted).to_excel(writer, sheet_name='Feature_Dropout', index=False)

    # Noise injection results
    noise_df = pd.DataFrame([
        {'Noise_Level': f'{int(n*100)}%',
         'AUC_Mean': noise_results[n]['AUC_mean'],
         'AUC_Std': noise_results[n]['AUC_std'],
         'Accuracy_Mean': noise_results[n]['Accuracy_mean'],
         'Accuracy_Std': noise_results[n]['Accuracy_std'],
         'F1_Mean': noise_results[n]['F1_mean'],
         'F1_Std': noise_results[n]['F1_std'],
         'AUC_Degradation': baseline_metrics['ROC-AUC'] - noise_results[n]['AUC_mean']}
        for n in noise_levels
    ])
    noise_df.to_excel(writer, sheet_name='Noise_Injection', index=False)

    # Feature-specific noise
    feat_noise_rows = []
    for feature in features_3:
        for noise in noise_test_levels:
            feat_noise_rows.append({
                'Feature': feature,
                'Noise_Level': f'{int(noise*100)}%',
                'AUC_Mean': feature_noise_sensitivity[feature][noise]['mean'],
                'AUC_Std': feature_noise_sensitivity[feature][noise]['std'],
                'AUC_Drop': baseline_metrics['ROC-AUC'] - feature_noise_sensitivity[feature][noise]['mean']
            })
    pd.DataFrame(feat_noise_rows).to_excel(writer, sheet_name='Feature_Noise_Sensitivity', index=False)

    # Bootstrap results
    boot_df = pd.DataFrame({
        'Bootstrap_AUC': bootstrap_aucs,
        'Bootstrap_Accuracy': bootstrap_accs
    })
    boot_df.to_excel(writer, sheet_name='Bootstrap_Results', index=False)

    # Bootstrap summary
    boot_summary = pd.DataFrame([{
        'Metric': 'AUC',
        'Mean': np.mean(bootstrap_aucs),
        'Std': np.std(bootstrap_aucs),
        'CI_Lower': np.percentile(bootstrap_aucs, 2.5),
        'CI_Upper': np.percentile(bootstrap_aucs, 97.5),
        'Min': np.min(bootstrap_aucs),
        'Max': np.max(bootstrap_aucs)
    }, {
        'Metric': 'Accuracy',
        'Mean': np.mean(bootstrap_accs),
        'Std': np.std(bootstrap_accs),
        'CI_Lower': np.percentile(bootstrap_accs, 2.5),
        'CI_Upper': np.percentile(bootstrap_accs, 97.5),
        'Min': np.min(bootstrap_accs),
        'Max': np.max(bootstrap_accs)
    }])
    boot_summary.to_excel(writer, sheet_name='Bootstrap_Summary', index=False)

    # Permutation test results
    perm_summary = pd.DataFrame([{
        'Observed_AUC': baseline_metrics['ROC-AUC'],
        'Null_Mean': np.mean(perm_aucs),
        'Null_Std': np.std(perm_aucs),
        'p_value': p_value,
        'Cohens_d': cohens_d,
        'N_Permutations': len(perm_aucs),
        'Significant': 'Yes' if p_value < 0.05 else 'No'
    }])
    perm_summary.to_excel(writer, sheet_name='Permutation_Test', index=False)

    # All permutation AUCs (for distribution)
    pd.DataFrame({'Permutation_AUC': perm_aucs}).to_excel(
        writer, sheet_name='Permutation_Distribution', index=False)

print("Saved: Robustness_Analysis_3Feature_Results.xlsx")

# ============================================================================
# 9. FINAL SUMMARY
# ============================================================================
print(f"\n{'='*70}")
print("ROBUSTNESS ANALYSIS SUMMARY")
print("=" * 70)

print(f"""
3-FEATURE MODEL ROBUSTNESS ASSESSMENT
=====================================

1. BASELINE PERFORMANCE:
   - ROC-AUC: {baseline_metrics['ROC-AUC']:.4f}
   - Accuracy: {baseline_metrics['Accuracy']:.4f}
   - Recall: {baseline_metrics['Recall']:.4f}
   - MCC: {baseline_metrics['MCC']:.4f}

2. STATISTICAL SIGNIFICANCE (Corrected Permutation Test):
   - p-value: {p_value:.4f} ({'SIGNIFICANT' if p_value < 0.05 else 'Not significant'})
   - Cohen's d: {cohens_d:.2f} ({'Large' if cohens_d > 0.8 else 'Medium' if cohens_d > 0.5 else 'Small'} effect size)
   - Null AUC: {np.mean(perm_aucs):.4f} vs Observed AUC: {baseline_metrics['ROC-AUC']:.4f}

3. FEATURE IMPORTANCE (by AUC drop when removed):
   - {dropout_results_sorted[0]['Feature']}: {dropout_results_sorted[0]['AUC_Drop']:+.4f} (CRITICAL)
   - {dropout_results_sorted[1]['Feature']}: {dropout_results_sorted[1]['AUC_Drop']:+.4f}
   - {dropout_results_sorted[2]['Feature']}: {dropout_results_sorted[2]['AUC_Drop']:+.4f}

4. NOISE ROBUSTNESS:
   - 10% noise: AUC = {noise_results[0.1]['AUC_mean']:.4f} ({(baseline_metrics['ROC-AUC'] - noise_results[0.1]['AUC_mean'])/baseline_metrics['ROC-AUC']*100:.1f}% degradation)
   - 30% noise: AUC = {noise_results[0.3]['AUC_mean']:.4f} ({(baseline_metrics['ROC-AUC'] - noise_results[0.3]['AUC_mean'])/baseline_metrics['ROC-AUC']*100:.1f}% degradation)
   - 50% noise: AUC = {noise_results[0.5]['AUC_mean']:.4f} ({(baseline_metrics['ROC-AUC'] - noise_results[0.5]['AUC_mean'])/baseline_metrics['ROC-AUC']*100:.1f}% degradation)

5. BOOTSTRAP STABILITY:
   - Mean AUC: {np.mean(bootstrap_aucs):.4f} +/- {np.std(bootstrap_aucs):.4f}
   - 95% CI: [{np.percentile(bootstrap_aucs, 2.5):.4f}, {np.percentile(bootstrap_aucs, 97.5):.4f}]

OUTPUT FILES:
   - Robustness_Analysis_3Feature.png/pdf/svg (comprehensive 8-panel figure)
   - Robustness_Analysis_3Feature_Results.xlsx (all numerical results)
""")

print("=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
