"""
Monte Carlo Robustness Analysis for 4-Feature PFS Classification Model
=======================================================================

Comprehensive robustness testing of the champion Naive Bayes classifier
for progression-free survival (PFS) prediction using four proteomic features.

Analyses performed:
    1. Targeted gene dropout (leave-one-out feature analysis)
    2. Random gene dropout at multiple dropout rates
    3. Expression noise injection (Gaussian noise scaled by feature std)
    4. Monte Carlo label shuffling (corrected permutation test)
    5. Feature-specific noise sensitivity profiling

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
                             precision_score, recall_score, confusion_matrix)
from sklearn.naive_bayes import GaussianNB
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import cross_val_score, StratifiedKFold
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
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
plt.style.use('seaborn-v0_8-whitegrid')

# TrueType font settings for Adobe Illustrator compatibility
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['axes.unicode_minus'] = False

# ============================================================================
# Load Data and Model
# ============================================================================
print("=" * 70)
print("MONTE CARLO ROBUSTNESS ANALYSIS FOR PFS CLASSIFICATION MODEL")
print("=" * 70)

with open(os.path.join(DATA_DIR, 'champion_model_NB_4features.pkl'), 'rb') as f:
    model_data = pickle.load(f)

model = model_data['model']
features = model_data['features']
print(f"\nModel: {model_data['algo']} with {len(features)} features")
print(f"Features: {features}")

train_df = pd.read_excel(os.path.join(DATA_DIR, 'PFS_ML_train_cof_wo_IG.xlsx'))
test_df = pd.read_excel(os.path.join(DATA_DIR, 'PFS_ML_test_cof_wo_IG.xlsx'))

# Encode PFS_group: L=0 (Long PFS), S=1 (Short PFS)
label_map = {'L': 0, 'S': 1}
train_df['PFS_group_encoded'] = train_df['PFS_group'].map(label_map)
test_df['PFS_group_encoded'] = test_df['PFS_group'].map(label_map)

X_train = train_df[features].values
y_train = train_df['PFS_group_encoded'].values
X_test = test_df[features].values
y_test = test_df['PFS_group_encoded'].values

print(f"\nTraining samples: {len(X_train)} ({sum(y_train)} Short PFS, {len(y_train)-sum(y_train)} Long PFS)")
print(f"Test samples: {len(X_test)} ({sum(y_test)} Short PFS, {len(y_test)-sum(y_test)} Long PFS)")

# Baseline performance
baseline_pred = model.predict(X_test)
baseline_proba = model.predict_proba(X_test)[:, 1]
baseline_metrics = {
    'Accuracy': accuracy_score(y_test, baseline_pred),
    'F1-Score': f1_score(y_test, baseline_pred),
    'ROC-AUC': roc_auc_score(y_test, baseline_proba),
    'Precision': precision_score(y_test, baseline_pred),
    'Recall': recall_score(y_test, baseline_pred)
}
print(f"\nBaseline Test Performance:")
for metric, value in baseline_metrics.items():
    print(f"  {metric}: {value:.4f}")

# ============================================================================
# 1. TARGETED GENE DROPOUT (Leave-One-Out Feature Analysis)
# ============================================================================
print("\n" + "=" * 70)
print("1. TARGETED GENE DROPOUT ANALYSIS")
print("=" * 70)
print("Assessing model dependency on each individual gene/protein")

targeted_dropout_results = []

for i, dropped_feature in enumerate(features):
    remaining_features = [f for f in features if f != dropped_feature]

    X_train_dropped = train_df[remaining_features].values
    X_test_dropped = test_df[remaining_features].values

    nb_dropped = GaussianNB()
    calibrated_dropped = CalibratedClassifierCV(nb_dropped, cv=3, method='sigmoid')
    calibrated_dropped.fit(X_train_dropped, y_train)

    pred_dropped = calibrated_dropped.predict(X_test_dropped)
    proba_dropped = calibrated_dropped.predict_proba(X_test_dropped)[:, 1]

    metrics = {
        'Dropped_Feature': dropped_feature,
        'Remaining_Features': remaining_features,
        'Accuracy': accuracy_score(y_test, pred_dropped),
        'F1-Score': f1_score(y_test, pred_dropped),
        'ROC-AUC': roc_auc_score(y_test, proba_dropped),
        'Accuracy_Drop': baseline_metrics['Accuracy'] - accuracy_score(y_test, pred_dropped),
        'F1_Drop': baseline_metrics['F1-Score'] - f1_score(y_test, pred_dropped),
        'AUC_Drop': baseline_metrics['ROC-AUC'] - roc_auc_score(y_test, proba_dropped)
    }
    targeted_dropout_results.append(metrics)

    print(f"\n  Dropped: {dropped_feature}")
    print(f"    Accuracy: {metrics['Accuracy']:.4f} (Delta = {metrics['Accuracy_Drop']:+.4f})")
    print(f"    F1-Score: {metrics['F1-Score']:.4f} (Delta = {metrics['F1_Drop']:+.4f})")
    print(f"    ROC-AUC:  {metrics['ROC-AUC']:.4f} (Delta = {metrics['AUC_Drop']:+.4f})")

targeted_df = pd.DataFrame(targeted_dropout_results)

# Feature importance ranking based on performance drop
print("\n  Feature Importance Ranking (by AUC drop when removed):")
ranked = targeted_df.sort_values('AUC_Drop', ascending=False)
for idx, row in ranked.iterrows():
    importance = row['AUC_Drop'] / baseline_metrics['ROC-AUC'] * 100
    print(f"    {row['Dropped_Feature']}: {row['AUC_Drop']:.4f} ({importance:.1f}% relative importance)")

# ============================================================================
# 2. RANDOM GENE DROPOUT ANALYSIS
# ============================================================================
print("\n" + "=" * 70)
print("2. RANDOM GENE DROPOUT ANALYSIS")
print("=" * 70)
print("Randomly dropping 1, 2, or 3 features to assess stability")

np.random.seed(42)
n_iterations = 100

random_dropout_results = {1: [], 2: [], 3: []}

for n_drop in [1, 2, 3]:
    print(f"\n  Dropping {n_drop} feature(s) randomly ({n_iterations} iterations)...")

    for iteration in range(n_iterations):
        drop_indices = np.random.choice(len(features), n_drop, replace=False)
        remaining_features = [f for i, f in enumerate(features) if i not in drop_indices]
        dropped_features = [features[i] for i in drop_indices]

        if len(remaining_features) == 0:
            continue

        X_train_dropped = train_df[remaining_features].values
        X_test_dropped = test_df[remaining_features].values

        nb_dropped = GaussianNB()
        if len(remaining_features) >= 2:
            calibrated_dropped = CalibratedClassifierCV(nb_dropped, cv=3, method='sigmoid')
        else:
            calibrated_dropped = CalibratedClassifierCV(nb_dropped, cv=2, method='sigmoid')

        try:
            calibrated_dropped.fit(X_train_dropped, y_train)
            pred_dropped = calibrated_dropped.predict(X_test_dropped)
            proba_dropped = calibrated_dropped.predict_proba(X_test_dropped)[:, 1]

            random_dropout_results[n_drop].append({
                'iteration': iteration,
                'dropped': dropped_features,
                'accuracy': accuracy_score(y_test, pred_dropped),
                'f1': f1_score(y_test, pred_dropped),
                'auc': roc_auc_score(y_test, proba_dropped)
            })
        except Exception:
            continue

    aucs = [r['auc'] for r in random_dropout_results[n_drop]]
    accs = [r['accuracy'] for r in random_dropout_results[n_drop]]

    print(f"    AUC: {np.mean(aucs):.4f} +/- {np.std(aucs):.4f} (range: {np.min(aucs):.4f} - {np.max(aucs):.4f})")
    print(f"    Accuracy: {np.mean(accs):.4f} +/- {np.std(accs):.4f}")

# ============================================================================
# 3. EXPRESSION NOISE INJECTION
# ============================================================================
print("\n" + "=" * 70)
print("3. EXPRESSION NOISE INJECTION ANALYSIS")
print("=" * 70)
print("Adding Gaussian noise to expression values to assess robustness")

noise_levels = [0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50]
n_noise_iterations = 100

noise_results = {level: {'accuracy': [], 'f1': [], 'auc': []} for level in noise_levels}

for noise_level in noise_levels:
    print(f"\n  Noise level: {noise_level*100:.0f}% of feature std...")

    for iteration in range(n_noise_iterations):
        # Noise scaled by per-feature standard deviation from training data
        feature_stds = np.std(X_train, axis=0)
        noise = np.random.normal(0, feature_stds * noise_level, X_test.shape)
        X_test_noisy = X_test + noise

        try:
            pred_noisy = model.predict(X_test_noisy)
            proba_noisy = model.predict_proba(X_test_noisy)[:, 1]

            noise_results[noise_level]['accuracy'].append(accuracy_score(y_test, pred_noisy))
            noise_results[noise_level]['f1'].append(f1_score(y_test, pred_noisy))
            noise_results[noise_level]['auc'].append(roc_auc_score(y_test, proba_noisy))
        except Exception:
            continue

    mean_auc = np.mean(noise_results[noise_level]['auc'])
    std_auc = np.std(noise_results[noise_level]['auc'])
    degradation = (baseline_metrics['ROC-AUC'] - mean_auc) / baseline_metrics['ROC-AUC'] * 100

    print(f"    AUC: {mean_auc:.4f} +/- {std_auc:.4f} ({degradation:.1f}% degradation from baseline)")

# ============================================================================
# 4. MONTE CARLO LABEL SHUFFLING (Permutation Test)
# ============================================================================
print("\n" + "=" * 70)
print("4. MONTE CARLO LABEL SHUFFLING (Corrected Permutation Test)")
print("=" * 70)
print("Shuffling TRAINING labels only, evaluating against TRUE test labels")
print("This tests whether the model learned real signal vs. noise")

n_permutations = 1000
permutation_aucs = []
permutation_accs = []
permutation_f1s = []

print(f"\n  Running {n_permutations} permutations...")

for i in range(n_permutations):
    if (i + 1) % 200 == 0:
        print(f"    Completed {i + 1}/{n_permutations} permutations...")

    # Shuffle ONLY training labels (correct methodology)
    y_train_shuffled = np.random.permutation(y_train)

    nb_perm = GaussianNB()
    calibrated_perm = CalibratedClassifierCV(nb_perm, cv=3, method='sigmoid')

    try:
        calibrated_perm.fit(X_train, y_train_shuffled)
        pred_perm = calibrated_perm.predict(X_test)
        proba_perm = calibrated_perm.predict_proba(X_test)[:, 1]

        # Evaluate against TRUE test labels (NOT shuffled)
        permutation_aucs.append(roc_auc_score(y_test, proba_perm))
        permutation_accs.append(accuracy_score(y_test, pred_perm))
        permutation_f1s.append(f1_score(y_test, pred_perm, zero_division=0))
    except Exception:
        continue

# Phipson-Smyth corrected empirical p-value: p = (r + 1) / (n + 1)
# Prevents p = 0.000 from finite permutation samples
empirical_p_auc = (np.sum(np.array(permutation_aucs) >= baseline_metrics['ROC-AUC']) + 1) / (len(permutation_aucs) + 1)
empirical_p_acc = (np.sum(np.array(permutation_accs) >= baseline_metrics['Accuracy']) + 1) / (len(permutation_accs) + 1)
empirical_p_f1 = (np.sum(np.array(permutation_f1s) >= baseline_metrics['F1-Score']) + 1) / (len(permutation_f1s) + 1)

print(f"\n  Null Distribution Statistics:")
print(f"    AUC:  mean={np.mean(permutation_aucs):.4f} +/- {np.std(permutation_aucs):.4f}")
print(f"    Accuracy: mean={np.mean(permutation_accs):.4f} +/- {np.std(permutation_accs):.4f}")
print(f"    F1-Score: mean={np.mean(permutation_f1s):.4f} +/- {np.std(permutation_f1s):.4f}")

print(f"\n  Observed vs Null:")
print(f"    AUC:  Observed={baseline_metrics['ROC-AUC']:.4f}, Null={np.mean(permutation_aucs):.4f}, p={empirical_p_auc:.4f}")
print(f"    Accuracy: Observed={baseline_metrics['Accuracy']:.4f}, Null={np.mean(permutation_accs):.4f}, p={empirical_p_acc:.4f}")
print(f"    F1-Score: Observed={baseline_metrics['F1-Score']:.4f}, Null={np.mean(permutation_f1s):.4f}, p={empirical_p_f1:.4f}")

# Effect size (Cohen's d)
cohens_d_auc = (baseline_metrics['ROC-AUC'] - np.mean(permutation_aucs)) / np.std(permutation_aucs)
print(f"\n  Effect Size (Cohen's d for AUC): {cohens_d_auc:.2f}")

# ============================================================================
# 5. COMPREHENSIVE VISUALIZATION
# ============================================================================
print("\n" + "=" * 70)
print("5. GENERATING COMPREHENSIVE VISUALIZATIONS")
print("=" * 70)

fig = plt.figure(figsize=(16, 12))

# Panel A: Targeted Gene Dropout
ax1 = fig.add_subplot(2, 3, 1)
x_pos = np.arange(len(features))
auc_drops = [targeted_df[targeted_df['Dropped_Feature'] == f]['AUC_Drop'].values[0] for f in features]
colors = ['#d62728' if d > 0.05 else '#2ca02c' if d < 0.02 else '#ff7f0e' for d in auc_drops]
bars = ax1.bar(x_pos, auc_drops, color=colors, edgecolor='black', linewidth=0.5)
ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
ax1.set_xticks(x_pos)
ax1.set_xticklabels(features, rotation=45, ha='right')
ax1.set_ylabel('AUC Drop (Delta)')
ax1.set_title('A. Targeted Gene Dropout\n(Feature Importance)')
ax1.set_ylim(-0.1, max(auc_drops) + 0.05)
for i, (bar, drop) in enumerate(zip(bars, auc_drops)):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
             f'{drop:.3f}', ha='center', va='bottom', fontsize=9)

# Panel B: Random Dropout Distribution
ax2 = fig.add_subplot(2, 3, 2)
dropout_data = []
labels = []
for n_drop in [1, 2, 3]:
    aucs = [r['auc'] for r in random_dropout_results[n_drop]]
    dropout_data.append(aucs)
    labels.append(f'{n_drop} feature(s)')

bp = ax2.boxplot(dropout_data, patch_artist=True, labels=labels)
colors_box = ['#1f77b4', '#ff7f0e', '#d62728']
for patch, color in zip(bp['boxes'], colors_box):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax2.axhline(y=baseline_metrics['ROC-AUC'], color='green', linestyle='--',
             linewidth=2, label=f'Baseline ({baseline_metrics["ROC-AUC"]:.3f})')
ax2.axhline(y=0.5, color='red', linestyle=':', linewidth=1.5, label='Random (0.5)')
ax2.set_ylabel('ROC-AUC')
ax2.set_xlabel('Number of Features Dropped')
ax2.set_title('B. Random Feature Dropout\n(Model Stability)')
ax2.legend(loc='lower left', fontsize=8)

# Panel C: Noise Injection
ax3 = fig.add_subplot(2, 3, 3)
noise_means = [np.mean(noise_results[level]['auc']) for level in noise_levels]
noise_stds = [np.std(noise_results[level]['auc']) for level in noise_levels]
noise_pct = [level * 100 for level in noise_levels]

ax3.errorbar(noise_pct, noise_means, yerr=noise_stds, fmt='o-', color='#1f77b4',
             capsize=4, capthick=1.5, markersize=8, linewidth=2)
ax3.axhline(y=baseline_metrics['ROC-AUC'], color='green', linestyle='--',
             linewidth=2, label=f'Baseline ({baseline_metrics["ROC-AUC"]:.3f})')
ax3.axhline(y=0.5, color='red', linestyle=':', linewidth=1.5, label='Random (0.5)')
ax3.fill_between(noise_pct,
                  [m - s for m, s in zip(noise_means, noise_stds)],
                  [m + s for m, s in zip(noise_means, noise_stds)],
                  alpha=0.2, color='#1f77b4')
ax3.set_xlabel('Noise Level (% of Feature Std)')
ax3.set_ylabel('ROC-AUC')
ax3.set_title('C. Expression Noise Injection\n(Measurement Error Tolerance)')
ax3.legend(loc='lower left', fontsize=8)
ax3.set_xlim(0, max(noise_pct) + 5)

# Panel D: Permutation Test Histogram
ax4 = fig.add_subplot(2, 3, 4)
ax4.hist(permutation_aucs, bins=40, color='lightblue', edgecolor='black',
         linewidth=0.5, alpha=0.7, label='Null Distribution')
ax4.axvline(x=baseline_metrics['ROC-AUC'], color='red', linestyle='-',
             linewidth=2.5, label=f'Observed ({baseline_metrics["ROC-AUC"]:.3f})')
ax4.axvline(x=np.mean(permutation_aucs), color='blue', linestyle='--',
             linewidth=1.5, label=f'Null Mean ({np.mean(permutation_aucs):.3f})')
ax4.axvline(x=0.5, color='gray', linestyle=':', linewidth=1.5, label='Random (0.5)')

p_text = f'p < 0.001' if empirical_p_auc < 0.001 else f'p = {empirical_p_auc:.3f}'
ax4.text(0.95, 0.95, p_text, transform=ax4.transAxes, fontsize=11,
         verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

ax4.set_xlabel('ROC-AUC')
ax4.set_ylabel('Frequency')
ax4.set_title('D. Monte Carlo Label Shuffling\n(Corrected Permutation Test)')
ax4.legend(loc='upper left', fontsize=8)

# Panel E: Summary Statistics Table
ax5 = fig.add_subplot(2, 3, 5)
ax5.axis('off')

summary_data = [
    ['Metric', 'Observed', 'Null Mean', 'Null Std', 'p-value', "Cohen's d"],
    ['ROC-AUC', f'{baseline_metrics["ROC-AUC"]:.4f}', f'{np.mean(permutation_aucs):.4f}',
     f'{np.std(permutation_aucs):.4f}', f'{empirical_p_auc:.4f}', f'{cohens_d_auc:.2f}'],
    ['Accuracy', f'{baseline_metrics["Accuracy"]:.4f}', f'{np.mean(permutation_accs):.4f}',
     f'{np.std(permutation_accs):.4f}', f'{empirical_p_acc:.4f}',
     f'{(baseline_metrics["Accuracy"] - np.mean(permutation_accs)) / np.std(permutation_accs):.2f}'],
    ['F1-Score', f'{baseline_metrics["F1-Score"]:.4f}', f'{np.mean(permutation_f1s):.4f}',
     f'{np.std(permutation_f1s):.4f}', f'{empirical_p_f1:.4f}',
     f'{(baseline_metrics["F1-Score"] - np.mean(permutation_f1s)) / np.std(permutation_f1s):.2f}']
]

table = ax5.table(cellText=summary_data, loc='center', cellLoc='center',
                  colWidths=[0.18, 0.15, 0.15, 0.15, 0.12, 0.15])
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.8)

for i in range(6):
    table[(0, i)].set_facecolor('#4472C4')
    table[(0, i)].set_text_props(color='white', weight='bold')

ax5.set_title('E. Statistical Summary\n(Observed vs. Null Distribution)', pad=20)

# Panel F: Feature Robustness Radar
ax6 = fig.add_subplot(2, 3, 6, polar=True)

robustness_metrics = []
for i, feature in enumerate(features):
    # Importance: normalised AUC drop when this feature is removed
    importance = auc_drops[i] / max(max(auc_drops), 0.01)

    # Noise tolerance: highest noise level at which mean AUC >= 0.8
    noise_tolerance = 0
    for level in noise_levels:
        if np.mean(noise_results[level]['auc']) >= 0.8:
            noise_tolerance = level
    noise_tolerance = noise_tolerance / max(noise_levels)

    robustness_metrics.append({
        'feature': feature,
        'importance': importance,
        'noise_tolerance': noise_tolerance
    })

categories = ['Feature\nImportance', 'Noise\nTolerance', 'Permutation\nSignificance',
              'Model\nStability']
N = len(categories)

noise_tol_avg = np.mean([r['noise_tolerance'] for r in robustness_metrics])
stability_score = np.mean([np.mean([r['auc'] for r in random_dropout_results[1]]),
                           np.mean([r['auc'] for r in random_dropout_results[2]])]) / baseline_metrics['ROC-AUC']

values = [
    np.mean(auc_drops) / baseline_metrics['ROC-AUC'],
    noise_tol_avg,
    1 - empirical_p_auc,
    stability_score
]
values += values[:1]  # Complete the loop

angles = [n / float(N) * 2 * np.pi for n in range(N)]
angles += angles[:1]

ax6.plot(angles, values, 'o-', linewidth=2, color='#2ca02c')
ax6.fill(angles, values, alpha=0.25, color='#2ca02c')
ax6.set_xticks(angles[:-1])
ax6.set_xticklabels(categories, size=9)
ax6.set_ylim(0, 1)
ax6.set_title('F. Model Robustness Profile', pad=15)

plt.suptitle('Monte Carlo Robustness Analysis\nPFS Classification Model (NB with 4 Features)',
             fontsize=14, fontweight='bold', y=1.02)

plt.tight_layout()

plt.savefig(os.path.join(OUTPUT_DIR, 'Monte_Carlo_Robustness_Analysis.png'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(os.path.join(OUTPUT_DIR, 'Monte_Carlo_Robustness_Analysis.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(os.path.join(OUTPUT_DIR, 'Monte_Carlo_Robustness_Analysis.svg'),
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print("\n  Saved: Monte_Carlo_Robustness_Analysis.png/pdf/svg")

# ============================================================================
# 6. FEATURE-SPECIFIC NOISE SENSITIVITY
# ============================================================================
print("\n" + "=" * 70)
print("6. FEATURE-SPECIFIC NOISE SENSITIVITY")
print("=" * 70)

feature_noise_sensitivity = {}

for feat_idx, feature in enumerate(features):
    print(f"\n  Analyzing {feature}...")

    sensitivities = []
    for noise_level in [0.1, 0.2, 0.3]:
        aucs_single_noise = []

        for _ in range(50):
            X_test_noisy = X_test.copy()
            feature_std = np.std(X_train[:, feat_idx])
            noise = np.random.normal(0, feature_std * noise_level, X_test.shape[0])
            X_test_noisy[:, feat_idx] += noise

            try:
                proba_noisy = model.predict_proba(X_test_noisy)[:, 1]
                aucs_single_noise.append(roc_auc_score(y_test, proba_noisy))
            except Exception:
                continue

        mean_auc = np.mean(aucs_single_noise)
        degradation = baseline_metrics['ROC-AUC'] - mean_auc
        sensitivities.append(degradation)

    feature_noise_sensitivity[feature] = {
        'mean_degradation': np.mean(sensitivities),
        'sensitivities': sensitivities
    }

    print(f"    Mean AUC degradation: {np.mean(sensitivities):.4f}")

fig2, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(features))
width = 0.25

for i, (level, label) in enumerate(zip([0.1, 0.2, 0.3], ['10% noise', '20% noise', '30% noise'])):
    degradations = [feature_noise_sensitivity[f]['sensitivities'][i] for f in features]
    bars = ax.bar(x + i*width, degradations, width, label=label, alpha=0.8)

ax.set_ylabel('AUC Degradation')
ax.set_xlabel('Feature')
ax.set_title('Feature-Specific Noise Sensitivity Analysis')
ax.set_xticks(x + width)
ax.set_xticklabels(features)
ax.legend()
ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)

plt.tight_layout()

plt.savefig(os.path.join(OUTPUT_DIR, 'Feature_Noise_Sensitivity.png'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(os.path.join(OUTPUT_DIR, 'Feature_Noise_Sensitivity.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(os.path.join(OUTPUT_DIR, 'Feature_Noise_Sensitivity.svg'),
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print("\n  Saved: Feature_Noise_Sensitivity.png/pdf/svg")

# ============================================================================
# 7. SAVE RESULTS TO EXCEL
# ============================================================================
print("\n" + "=" * 70)
print("7. SAVING RESULTS TO EXCEL")
print("=" * 70)

with pd.ExcelWriter(os.path.join(OUTPUT_DIR, 'Monte_Carlo_Robustness_Results.xlsx'),
                     engine='openpyxl') as writer:

    # Baseline metrics
    baseline_df = pd.DataFrame([baseline_metrics])
    baseline_df.to_excel(writer, sheet_name='Baseline_Metrics', index=False)

    # Targeted dropout
    targeted_df.to_excel(writer, sheet_name='Targeted_Dropout', index=False)

    # Random dropout summary
    random_dropout_summary = []
    for n_drop in [1, 2, 3]:
        aucs = [r['auc'] for r in random_dropout_results[n_drop]]
        random_dropout_summary.append({
            'Features_Dropped': n_drop,
            'Mean_AUC': np.mean(aucs),
            'Std_AUC': np.std(aucs),
            'Min_AUC': np.min(aucs),
            'Max_AUC': np.max(aucs),
            'N_Iterations': len(aucs)
        })
    pd.DataFrame(random_dropout_summary).to_excel(writer, sheet_name='Random_Dropout', index=False)

    # Noise injection
    noise_summary = []
    for level in noise_levels:
        noise_summary.append({
            'Noise_Level_Pct': level * 100,
            'Mean_AUC': np.mean(noise_results[level]['auc']),
            'Std_AUC': np.std(noise_results[level]['auc']),
            'Mean_Accuracy': np.mean(noise_results[level]['accuracy']),
            'Mean_F1': np.mean(noise_results[level]['f1']),
            'AUC_Degradation_Pct': (baseline_metrics['ROC-AUC'] - np.mean(noise_results[level]['auc'])) / baseline_metrics['ROC-AUC'] * 100
        })
    pd.DataFrame(noise_summary).to_excel(writer, sheet_name='Noise_Injection', index=False)

    # Permutation test
    perm_df = pd.DataFrame({
        'Iteration': range(len(permutation_aucs)),
        'AUC': permutation_aucs,
        'Accuracy': permutation_accs,
        'F1': permutation_f1s
    })
    perm_df.to_excel(writer, sheet_name='Permutation_Test', index=False)

    # Permutation summary
    perm_summary = pd.DataFrame([{
        'Metric': 'ROC-AUC',
        'Observed': baseline_metrics['ROC-AUC'],
        'Null_Mean': np.mean(permutation_aucs),
        'Null_Std': np.std(permutation_aucs),
        'Empirical_P_Value': empirical_p_auc,
        'Cohens_d': cohens_d_auc,
        'N_Permutations': len(permutation_aucs)
    }, {
        'Metric': 'Accuracy',
        'Observed': baseline_metrics['Accuracy'],
        'Null_Mean': np.mean(permutation_accs),
        'Null_Std': np.std(permutation_accs),
        'Empirical_P_Value': empirical_p_acc,
        'Cohens_d': (baseline_metrics['Accuracy'] - np.mean(permutation_accs)) / np.std(permutation_accs),
        'N_Permutations': len(permutation_accs)
    }, {
        'Metric': 'F1-Score',
        'Observed': baseline_metrics['F1-Score'],
        'Null_Mean': np.mean(permutation_f1s),
        'Null_Std': np.std(permutation_f1s),
        'Empirical_P_Value': empirical_p_f1,
        'Cohens_d': (baseline_metrics['F1-Score'] - np.mean(permutation_f1s)) / np.std(permutation_f1s),
        'N_Permutations': len(permutation_f1s)
    }])
    perm_summary.to_excel(writer, sheet_name='Permutation_Summary', index=False)

    # Feature noise sensitivity
    noise_sens_df = pd.DataFrame([
        {'Feature': f,
         '10%_Noise_Degradation': feature_noise_sensitivity[f]['sensitivities'][0],
         '20%_Noise_Degradation': feature_noise_sensitivity[f]['sensitivities'][1],
         '30%_Noise_Degradation': feature_noise_sensitivity[f]['sensitivities'][2],
         'Mean_Degradation': feature_noise_sensitivity[f]['mean_degradation']}
        for f in features
    ])
    noise_sens_df.to_excel(writer, sheet_name='Feature_Noise_Sensitivity', index=False)

print("  Saved: Monte_Carlo_Robustness_Results.xlsx")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY: MODEL ROBUSTNESS ASSESSMENT")
print("=" * 70)

print(f"""
MODEL: Naive Bayes with 4 features ({', '.join(features)})
BASELINE PERFORMANCE: AUC={baseline_metrics['ROC-AUC']:.4f}, Accuracy={baseline_metrics['Accuracy']:.4f}

1. TARGETED GENE DROPOUT (Feature Importance):
   - Most important feature: {ranked.iloc[0]['Dropped_Feature']} (AUC drop: {ranked.iloc[0]['AUC_Drop']:.4f})
   - Least important feature: {ranked.iloc[-1]['Dropped_Feature']} (AUC drop: {ranked.iloc[-1]['AUC_Drop']:.4f})

2. RANDOM GENE DROPOUT (Model Stability):
   - Dropping 1 feature: AUC = {np.mean([r['auc'] for r in random_dropout_results[1]]):.4f} +/- {np.std([r['auc'] for r in random_dropout_results[1]]):.4f}
   - Dropping 2 features: AUC = {np.mean([r['auc'] for r in random_dropout_results[2]]):.4f} +/- {np.std([r['auc'] for r in random_dropout_results[2]]):.4f}

3. NOISE TOLERANCE (Measurement Error Robustness):
   - At 10% noise: AUC = {np.mean(noise_results[0.10]['auc']):.4f} ({(baseline_metrics['ROC-AUC'] - np.mean(noise_results[0.10]['auc']))/baseline_metrics['ROC-AUC']*100:.1f}% degradation)
   - At 20% noise: AUC = {np.mean(noise_results[0.20]['auc']):.4f} ({(baseline_metrics['ROC-AUC'] - np.mean(noise_results[0.20]['auc']))/baseline_metrics['ROC-AUC']*100:.1f}% degradation)
   - At 30% noise: AUC = {np.mean(noise_results[0.30]['auc']):.4f} ({(baseline_metrics['ROC-AUC'] - np.mean(noise_results[0.30]['auc']))/baseline_metrics['ROC-AUC']*100:.1f}% degradation)

4. CORRECTED PERMUTATION TEST (Statistical Significance):
   - Methodology: Shuffle TRAINING labels only, evaluate against TRUE test labels
   - Observed AUC: {baseline_metrics['ROC-AUC']:.4f}
   - Null distribution: mean={np.mean(permutation_aucs):.4f} +/- {np.std(permutation_aucs):.4f}
   - p-value: {empirical_p_auc:.4f} ({'SIGNIFICANT' if empirical_p_auc < 0.05 else 'Not significant'} at alpha=0.05)
   - Effect size (Cohen's d): {cohens_d_auc:.2f} ({'Large' if abs(cohens_d_auc) > 0.8 else 'Medium' if abs(cohens_d_auc) > 0.5 else 'Small'} effect)

OUTPUT FILES:
   - Monte_Carlo_Robustness_Analysis.png/pdf/svg (comprehensive 6-panel figure)
   - Feature_Noise_Sensitivity.png/pdf/svg (feature-specific analysis)
   - Monte_Carlo_Robustness_Results.xlsx (all numerical results)
""")

print("=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
