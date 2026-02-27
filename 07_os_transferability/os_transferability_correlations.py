# -*- coding: utf-8 -*-
"""
OS Transferability - Correlation Analysis
==========================================
Tests whether PFS-trained proteomic signatures (NB, Cox, ensemble) correlate
with Overall Survival, providing cross-endpoint transferability evidence.

Analyses:
1. Per-feature Spearman correlations with OS (FDR-corrected)
2. PFS vs OS correlation (patient-level)
3. Cox risk score, NB probability, and ensemble probability vs OS
4. PFS-OS group concordance (Cohen's kappa)
5. Feature-level boxplots: Short OS vs Long OS (t-test, nominal)

Reference:
Khoshnevis et al. Molecular Cancer (2026).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.naive_bayes import GaussianNB
from sklearn.calibration import CalibratedClassifierCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import cohen_kappa_score
import warnings
import os

warnings.filterwarnings('ignore')

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# PUBLICATION-QUALITY FIGURE SETTINGS
# ============================================================================
plt.rcParams['pdf.fonttype'] = 42   # TrueType (editable in Illustrator)
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

DPI = 300
FORMATS = ['png', 'pdf']

# ============================================================================
# MODEL PARAMETERS (PFS-trained, transferred to OS)
# ============================================================================
NB_FEATURES = ['DYNC2H1', 'ECM2', 'PPIB']
COX_FEATURES = ['TFRC', 'APOF', 'ANG', 'FABP4']
ALL_FEATURES = NB_FEATURES + COX_FEATURES
COX_BETAS = {'TFRC': 0.390057, 'APOF': 0.391837, 'ANG': 0.269207, 'FABP4': 0.280084}

# Colors (matched to PFS ensemble figures)
COL_NB = '#2E86AB'
COL_COX = '#A23B72'
COL_ENS = '#28A745'
COL_SHORT = '#E74C3C'
COL_LONG = '#2E86AB'


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def savefig(name, fig=None):
    """Save figure in PNG and PDF formats."""
    if fig is None:
        fig = plt.gcf()
    for fmt in FORMATS:
        fig.savefig(os.path.join(OUTPUT_DIR, f"{name}.{fmt}"),
                    bbox_inches='tight', dpi=DPI, facecolor='white')
    plt.close(fig)
    print(f"  Saved: {name}")


def add_pval_text(ax, r, p, x=0.05, y=0.95, label='nominal'):
    """Add correlation statistics annotation to axis."""
    ptxt = 'p < 0.001' if p < 0.001 else f'p = {p:.3f}'
    ax.text(x, y, f'r = {r:.3f}\n{ptxt} ({label})', transform=ax.transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      alpha=0.8, edgecolor='gray'))


# ============================================================================
# LOAD DATA
# ============================================================================
print("=" * 70)
print("OS TRANSFERABILITY - CORRELATION ANALYSIS")
print("=" * 70)

# OS-labelled dataset (33 patients with available OS data)
df = pd.read_excel(os.path.join(DATA_DIR, 'FinalwinnerML_Cox_dataset_OS.xlsx'))

# Full PFS dataset (for NB training on the original PFS endpoint)
df_pfs = pd.read_excel(os.path.join(DATA_DIR, 'FinalwinnerML_Cox_dataset.xlsx'))

# Merge PFS information
df = df.merge(df_pfs[['patient_ID', 'PFS', 'PFS_group']], on='patient_ID', how='left')

print(f"Patients: {len(df)}")
print(f"OS_group: S={sum(df['OS_group'] == 'S')}, L={sum(df['OS_group'] == 'L')}")
print(f"OS range: {df['OS'].min():.0f} - {df['OS'].max():.0f} days")

# ============================================================================
# COMPUTE MODEL SCORES
# ============================================================================
# Cox risk score (PFS-trained betas applied as fixed coefficients)
df['Cox_risk'] = sum(COX_BETAS[f] * df[f] for f in COX_FEATURES)
cox_median = df['Cox_risk'].median()
df['Cox_prob'] = 1 / (1 + np.exp(-(df['Cox_risk'] - cox_median)))

# NB probability (trained on full PFS cohort, applied to OS subset)
scaler = StandardScaler()
X_nb_scaled = scaler.fit_transform(df_pfs[NB_FEATURES].values)
y_pfs = df_pfs['PFS_group'].map({'L': 0, 'S': 1}).values
base_nb = GaussianNB()
cal_nb = CalibratedClassifierCV(base_nb, method='sigmoid', cv=3)
cal_nb.fit(X_nb_scaled, y_pfs)

X_nb_os = scaler.transform(df[NB_FEATURES].values)
df['NB_prob'] = cal_nb.predict_proba(X_nb_os)[:, 1]

# Ensemble probability (weighted_1_2: Cox-weighted)
df['Ens_prob'] = (df['NB_prob'] + 2 * df['Cox_prob']) / 3

# Binary labels
df['y_os'] = df['OS_group'].map({'L': 0, 'S': 1})

print(f"\nCox risk score range: {df['Cox_risk'].min():.2f} - {df['Cox_risk'].max():.2f}")
print(f"NB prob range: {df['NB_prob'].min():.3f} - {df['NB_prob'].max():.3f}")
print(f"Ensemble prob range: {df['Ens_prob'].min():.3f} - {df['Ens_prob'].max():.3f}")

# ============================================================================
# FIGURE: Multi-panel correlation overview (3x3 grid, Fig S4)
# ============================================================================
print("\n[1] Multi-panel correlation figure...")

fig = plt.figure(figsize=(16, 12))
gs = GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.35)

# --- Panel A: PFS vs OS scatter ---
ax_a = fig.add_subplot(gs[0, 0])
colors_a = [COL_SHORT if g == 'S' else COL_LONG for g in df['OS_group']]
ax_a.scatter(df['PFS'], df['OS'], c=colors_a, s=50, alpha=0.7,
             edgecolors='black', linewidth=0.5, zorder=3)
r_pfs_os, p_pfs_os = stats.spearmanr(df['PFS'], df['OS'])
z = np.polyfit(df['PFS'], df['OS'], 1)
x_line = np.linspace(df['PFS'].min(), df['PFS'].max(), 100)
ax_a.plot(x_line, np.polyval(z, x_line), '--', color='gray', alpha=0.7, linewidth=1.5)
add_pval_text(ax_a, r_pfs_os, p_pfs_os)
ax_a.set_xlabel('PFS (days)')
ax_a.set_ylabel('OS (days)')
ax_a.set_title('A. PFS vs OS Correlation', fontweight='bold', fontsize=11)
ax_a.legend(handles=[mpatches.Patch(color=COL_SHORT, label='Short OS'),
                      mpatches.Patch(color=COL_LONG, label='Long OS')],
            fontsize=8, loc='lower right')

# --- Panel B: Cox risk score vs OS ---
ax_b = fig.add_subplot(gs[0, 1])
ax_b.scatter(df['Cox_risk'], df['OS'], c=colors_a, s=50, alpha=0.7,
             edgecolors='black', linewidth=0.5, zorder=3)
r_cox_os, p_cox_os = stats.spearmanr(df['Cox_risk'], df['OS'])
z2 = np.polyfit(df['Cox_risk'], df['OS'], 1)
x_line2 = np.linspace(df['Cox_risk'].min(), df['Cox_risk'].max(), 100)
ax_b.plot(x_line2, np.polyval(z2, x_line2), '--', color='gray', alpha=0.7, linewidth=1.5)
add_pval_text(ax_b, r_cox_os, p_cox_os)
ax_b.set_xlabel('Cox Risk Score (PFS-trained)')
ax_b.set_ylabel('OS (days)')
ax_b.set_title('B. Cox Risk Score vs OS', fontweight='bold', fontsize=11)

# --- Panel C: NB probability vs OS ---
ax_c = fig.add_subplot(gs[0, 2])
ax_c.scatter(df['NB_prob'], df['OS'], c=colors_a, s=50, alpha=0.7,
             edgecolors='black', linewidth=0.5, zorder=3)
r_nb_os, p_nb_os = stats.spearmanr(df['NB_prob'], df['OS'])
z3 = np.polyfit(df['NB_prob'], df['OS'], 1)
x_line3 = np.linspace(df['NB_prob'].min(), df['NB_prob'].max(), 100)
ax_c.plot(x_line3, np.polyval(z3, x_line3), '--', color='gray', alpha=0.7, linewidth=1.5)
add_pval_text(ax_c, r_nb_os, p_nb_os)
ax_c.set_xlabel('NB P(Short) (PFS-trained)')
ax_c.set_ylabel('OS (days)')
ax_c.set_title('C. NB Probability vs OS', fontweight='bold', fontsize=11)

# --- Panel D: Ensemble probability vs OS ---
ax_d = fig.add_subplot(gs[1, 0])
ax_d.scatter(df['Ens_prob'], df['OS'], c=colors_a, s=50, alpha=0.7,
             edgecolors='black', linewidth=0.5, zorder=3)
r_ens_os, p_ens_os = stats.spearmanr(df['Ens_prob'], df['OS'])
z4 = np.polyfit(df['Ens_prob'], df['OS'], 1)
x_line4 = np.linspace(df['Ens_prob'].min(), df['Ens_prob'].max(), 100)
ax_d.plot(x_line4, np.polyval(z4, x_line4), '--', color='gray', alpha=0.7, linewidth=1.5)
add_pval_text(ax_d, r_ens_os, p_ens_os)
ax_d.set_xlabel('Ensemble P(Short) (PFS-trained)')
ax_d.set_ylabel('OS (days)')
ax_d.set_title('D. Ensemble Probability vs OS', fontweight='bold', fontsize=11)

# --- Panel E: PFS-OS group concordance heatmap ---
ax_e = fig.add_subplot(gs[1, 1])
ct = pd.crosstab(df['PFS_group'], df['OS_group'])
ct_pct = ct / ct.sum().sum() * 100
ct = ct.reindex(index=['L', 'S'], columns=['L', 'S'])
ct_pct = ct_pct.reindex(index=['L', 'S'], columns=['L', 'S'])

im = ax_e.imshow(ct.values, cmap='Blues', aspect='auto', vmin=0, vmax=ct.values.max())
for i in range(2):
    for j in range(2):
        val = ct.values[i, j]
        pct = ct_pct.values[i, j]
        color = 'white' if val > ct.values.max() * 0.6 else 'black'
        ax_e.text(j, i, f'{val}\n({pct:.0f}%)', ha='center', va='center',
                  fontsize=11, fontweight='bold', color=color)
ax_e.set_xticks([0, 1])
ax_e.set_xticklabels(['Long OS', 'Short OS'])
ax_e.set_yticks([0, 1])
ax_e.set_yticklabels(['Long PFS', 'Short PFS'])
ax_e.set_xlabel('OS Group')
ax_e.set_ylabel('PFS Group')
kappa = cohen_kappa_score(df['PFS_group'], df['OS_group'])
concordance = sum(df['PFS_group'] == df['OS_group']) / len(df) * 100
ax_e.set_title(f'E. PFS-OS Concordance\n({concordance:.0f}%, \u03ba={kappa:.2f})',
               fontweight='bold', fontsize=11)

# --- Panel F: Per-feature correlations with OS (bar chart) ---
ax_f = fig.add_subplot(gs[1, 2])
feature_corrs = []
for feat in ALL_FEATURES:
    r, p = stats.spearmanr(df[feat], df['OS'])
    feature_corrs.append({'Feature': feat, 'r': r, 'p': p,
                           'Component': 'NB' if feat in NB_FEATURES else 'Cox'})
fc_df = pd.DataFrame(feature_corrs)

# FDR correction (Benjamini-Hochberg) for 7 per-feature tests
_, fc_df['p_fdr'], _, _ = multipletests(fc_df['p'].values, method='fdr_bh')

fc_df = fc_df.sort_values('r')

colors_bar = [COL_NB if c == 'NB' else COL_COX for c in fc_df['Component']]
ax_f.barh(range(len(fc_df)), fc_df['r'], color=colors_bar,
          edgecolor='black', linewidth=0.5, height=0.6)
ax_f.set_yticks(range(len(fc_df)))
ax_f.set_yticklabels(fc_df['Feature'], fontsize=9)
ax_f.set_xlabel('Spearman r with OS')
ax_f.axvline(x=0, color='black', linewidth=0.8, linestyle='-')
for i, (_, row) in enumerate(fc_df.iterrows()):
    q = row['p_fdr']
    star = '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'ns'
    if row['r'] >= 0:
        x_pos = row['r'] + 0.015
        ha = 'left'
    else:
        x_pos = min(row['r'] - 0.015, -0.06)
        ha = 'right'
    ax_f.text(x_pos, i, star, va='center', ha=ha, fontsize=8, color='black')
ax_f.legend(handles=[mpatches.Patch(color=COL_NB, label='NB features'),
                      mpatches.Patch(color=COL_COX, label='Cox features')],
            fontsize=8, loc='lower right')
ax_f.set_title('F. Feature-OS Correlations (FDR)', fontweight='bold', fontsize=11)

# --- Panels G-I: Feature boxplots for top 3 correlated ---
fc_df_sorted = fc_df.sort_values('p', ascending=True)
top3 = fc_df_sorted.head(3)

for idx, (panel_label, (_, row)) in enumerate(zip(['G', 'H', 'I'], top3.iterrows())):
    ax = fig.add_subplot(gs[2, idx])
    feat = row['Feature']
    short_vals = df.loc[df['OS_group'] == 'S', feat]
    long_vals = df.loc[df['OS_group'] == 'L', feat]

    bp = ax.boxplot([long_vals, short_vals], positions=[1, 2], widths=0.5,
                    patch_artist=True, showfliers=True)
    bp['boxes'][0].set_facecolor(COL_LONG)
    bp['boxes'][0].set_alpha(0.5)
    bp['boxes'][1].set_facecolor(COL_SHORT)
    bp['boxes'][1].set_alpha(0.5)
    for element in ['whiskers', 'caps', 'medians']:
        for item in bp[element]:
            item.set_color('black')

    # Individual points with jitter
    np.random.seed(idx)
    for val, pos in [(long_vals, 1), (short_vals, 2)]:
        jitter = np.random.normal(0, 0.04, size=len(val))
        color_pt = COL_LONG if pos == 1 else COL_SHORT
        ax.scatter(np.full(len(val), pos) + jitter, val, c=color_pt, s=20,
                   alpha=0.6, edgecolors='black', linewidth=0.3, zorder=3)

    # T-test (nominal, per-feature)
    t_stat, t_p = stats.ttest_ind(short_vals, long_vals)
    t_p_val = float(t_p)
    star = '***' if t_p_val < 0.001 else '**' if t_p_val < 0.01 else '*' if t_p_val < 0.05 else 'ns'
    y_max = max(float(short_vals.max()), float(long_vals.max()))
    y_min = min(float(short_vals.min()), float(long_vals.min()))
    y_range = y_max - y_min
    bracket_y = y_max + 0.05 * y_range
    ax.plot([1, 1, 2, 2],
            [bracket_y, bracket_y + 0.02 * y_range,
             bracket_y + 0.02 * y_range, bracket_y],
            color='black', linewidth=1)
    p_text = f'p={t_p_val:.3f} (nominal)' if t_p_val >= 0.001 else 'p<0.001 (nominal)'
    ax.text(1.5, bracket_y + 0.03 * y_range, f'{star}  {p_text}',
            ha='center', va='bottom', fontsize=8)
    ax.set_ylim(y_min - 0.05 * y_range, bracket_y + 0.20 * y_range)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Long OS', 'Short OS'])
    ax.set_ylabel(str(feat))
    comp = 'NB' if feat in NB_FEATURES else 'Cox'
    ax.set_title(f'{panel_label}. {feat} ({comp})\nr={row["r"]:.3f}',
                 fontweight='bold', fontsize=11)

fig.suptitle('OS Transferability: Correlation Evidence for PFS-Trained Proteomic Signatures',
             fontsize=14, fontweight='bold', y=1.01)
fig.tight_layout(rect=(0, 0, 1, 0.96))
savefig('OS_Transferability_Correlations', fig)

# ============================================================================
# SUMMARY STATISTICS TABLE
# ============================================================================
print("\n[2] Summary statistics...")

summary = []
feat_pvals = []
for feat in ALL_FEATURES:
    r, p = stats.spearmanr(df[feat], df['OS'])
    t_stat, t_p = stats.ttest_ind(
        df.loc[df['OS_group'] == 'S', feat],
        df.loc[df['OS_group'] == 'L', feat]
    )
    feat_pvals.append(float(p))
    summary.append({
        'Variable': feat, 'Type': 'NB Feature' if feat in NB_FEATURES else 'Cox Feature',
        'Spearman_r_OS': r, 'Spearman_p_nominal': float(p),
        'Spearman_p_FDR': np.nan,
        'GroupTest_p_nominal': float(t_p),
        'Test_type': 't-test',
        'Mean_ShortOS': df.loc[df['OS_group'] == 'S', feat].mean(),
        'Mean_LongOS': df.loc[df['OS_group'] == 'L', feat].mean()
    })

# FDR correction for per-feature Spearman p-values
_, fdr_qvals, _, _ = multipletests(feat_pvals, method='fdr_bh')
for i in range(len(ALL_FEATURES)):
    summary[i]['Spearman_p_FDR'] = fdr_qvals[i]

for score_col, label in [('Cox_risk', 'Cox Risk Score'), ('Cox_prob', 'Cox P(Short)'),
                          ('NB_prob', 'NB P(Short)'), ('Ens_prob', 'Ensemble P(Short)')]:
    r, p = stats.spearmanr(df[score_col], df['OS'])
    u, u_p = stats.mannwhitneyu(
        df.loc[df['OS_group'] == 'S', score_col],
        df.loc[df['OS_group'] == 'L', score_col], alternative='two-sided'
    )
    summary.append({
        'Variable': label, 'Type': 'Model Score',
        'Spearman_r_OS': r, 'Spearman_p_nominal': float(p),
        'Spearman_p_FDR': np.nan,
        'GroupTest_p_nominal': float(u_p),
        'Test_type': 'Mann-Whitney U',
        'Mean_ShortOS': df.loc[df['OS_group'] == 'S', score_col].mean(),
        'Mean_LongOS': df.loc[df['OS_group'] == 'L', score_col].mean()
    })

r_pfs, p_pfs = stats.spearmanr(df['PFS'], df['OS'])
summary.append({
    'Variable': 'PFS vs OS', 'Type': 'Endpoint Correlation',
    'Spearman_r_OS': r_pfs, 'Spearman_p_nominal': float(p_pfs),
    'Spearman_p_FDR': np.nan,
    'GroupTest_p_nominal': np.nan,
    'Test_type': 'N/A',
    'Mean_ShortOS': df.loc[df['OS_group'] == 'S', 'PFS'].mean(),
    'Mean_LongOS': df.loc[df['OS_group'] == 'L', 'PFS'].mean()
})

summary_df = pd.DataFrame(summary)

concordance_info = pd.DataFrame({
    'Metric': ['PFS-OS Group Concordance', 'Cohen Kappa', 'N patients',
               'Short OS (S)', 'Long OS (L)', 'OS median (days)'],
    'Value': [f'{concordance:.1f}%', f'{kappa:.3f}', len(df),
              sum(df['OS_group'] == 'S'), sum(df['OS_group'] == 'L'),
              f"{df['OS'].median():.0f}"]
})

outpath = os.path.join(OUTPUT_DIR, 'OS_Transferability_Correlation_Results.xlsx')
with pd.ExcelWriter(outpath, engine='openpyxl') as writer:
    summary_df.to_excel(writer, sheet_name='Correlation_Summary', index=False)
    concordance_info.to_excel(writer, sheet_name='Concordance_Info', index=False)

print(f"\nSaved: {outpath}")

# ============================================================================
# KEY RESULTS SUMMARY
# ============================================================================
print("\n" + "=" * 70)
print("KEY CORRELATION RESULTS")
print("=" * 70)
print(f"\nPFS vs OS: r={r_pfs_os:.3f}, p={p_pfs_os:.4f}")
print(f"Cox risk vs OS: r={r_cox_os:.3f}, p={p_cox_os:.4f}")
print(f"NB prob vs OS: r={r_nb_os:.3f}, p={p_nb_os:.4f}")
print(f"Ensemble prob vs OS: r={r_ens_os:.3f}, p={p_ens_os:.4f}")
print(f"\nPFS-OS concordance: {concordance:.1f}% (kappa={kappa:.3f})")
print(f"\nPer-feature correlations with OS (FDR-corrected):")
for _, row in fc_df.sort_values('p').iterrows():
    q = row['p_fdr']
    star = '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else 'ns'
    print(f"  {row['Feature']:>10s} ({row['Component']}): "
          f"r={row['r']:.3f}, p_nom={row['p']:.4f}, q_FDR={q:.4f} {star}")

print("\nDone!")
