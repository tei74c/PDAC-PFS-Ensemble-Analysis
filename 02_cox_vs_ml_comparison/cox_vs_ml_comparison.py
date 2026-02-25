# -*- coding: utf-8 -*-
"""
Cox vs ML Head-to-Head Comparison
==================================
Bootstrap cross-validation comparison of Cox PH and Naive Bayes classifiers.
Generates Figures 4A-G in the manuscript.

Reference:
Khoshnevis et al. Molecular Cancer (2026).
"""

# =============================================================================
# CONFIGURATION
# =============================================================================

import os

# ---------------------------------------------------------------------------
# FILE PATHS (relative to script location)
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

FULL_DATASET = os.path.join(DATA_DIR, "FinalwinnerML_Cox_dataset.xlsx")
NB_PICKLE = os.path.join(DATA_DIR, "champion_model_NB_3features.pkl")

# ---------------------------------------------------------------------------
# COX MODEL PARAMETERS (from Cox model_winner.xlsx)
# ---------------------------------------------------------------------------
COX_FEATURES = ['TFRC', 'APOF', 'ANG', 'FABP4']
COX_BETAS = {
    'TFRC': 0.390057,
    'APOF': 0.391837,
    'ANG': 0.269207,
    'FABP4': 0.280084
}

# ---------------------------------------------------------------------------
# ML MODEL PARAMETERS (3-Feature Model)
# ---------------------------------------------------------------------------
NB_FEATURES = ['DYNC2H1', 'ECM2', 'PPIB']

# ---------------------------------------------------------------------------
# CLASS CONFIGURATION
# ---------------------------------------------------------------------------
LABEL_COLUMN = 'PFS_group'
CLASS_MAP = {'L': 0, 'S': 1}  # L=Long PFS (negative), S=Short PFS (positive)
CLASS_NAMES = ['Long PFS (L)', 'Short PFS (S)']
POSITIVE_CLASS = 'S'
POSITIVE_LABEL = 1

# Survival columns (for Kaplan-Meier validation)
TIME_COLUMN = 'PFS'
EVENT_COLUMN = 'Event'

# ---------------------------------------------------------------------------
# BOOTSTRAP CV PARAMETERS
# ---------------------------------------------------------------------------
N_BOOTSTRAP = 200
CV_N_REPEATS = 20
CV_N_FOLDS = 5
TOTAL_EVALUATIONS = N_BOOTSTRAP * CV_N_REPEATS * CV_N_FOLDS  # 20,000

# ---------------------------------------------------------------------------
# STATISTICAL PARAMETERS
# ---------------------------------------------------------------------------
CI_LEVEL = 0.95
FDR_ALPHA = 0.05
N_PERMUTATIONS = 1000
DELONG_ALPHA = 0.05

# ---------------------------------------------------------------------------
# RISK STRATIFICATION
# ---------------------------------------------------------------------------
RISK_THRESHOLD_METHOD = 'median'
RISK_THRESHOLD_FIXED = 0.5

# ---------------------------------------------------------------------------
# DECISION CURVE ANALYSIS
# ---------------------------------------------------------------------------
DCA_THRESHOLD_RANGE = (0.01, 0.99)
DCA_N_THRESHOLDS = 100

# ---------------------------------------------------------------------------
# NRI/IDI PARAMETERS
# ---------------------------------------------------------------------------
NRI_THRESHOLD = 0.5

# ---------------------------------------------------------------------------
# VISUALIZATION PARAMETERS
# ---------------------------------------------------------------------------
FIGURE_DPI = 300
SAVE_FORMATS = ['png', 'pdf', 'svg']
COLORMAP = {
    'NB': '#2E86AB',
    'Cox': '#A23B72',
    'ZeroR': '#F18F01',
    'Champion': '#2E86AB',
    'Majority': 'lightcoral'
}

FONT_SIZE = 11
TITLE_SIZE = 14
LABEL_SIZE = 12
TICK_SIZE = 10

TRUETYPE_SETTINGS = {
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
}

# ---------------------------------------------------------------------------
# COMPUTATIONAL PARAMETERS
# ---------------------------------------------------------------------------
RANDOM_STATE = 42
N_JOBS = -1
VERBOSE = True

# ---------------------------------------------------------------------------
# METRICS TO COMPUTE
# ---------------------------------------------------------------------------
CORE_METRICS = [
    'Accuracy', 'Balanced_Accuracy', 'F1', 'Precision', 'Recall',
    'Specificity', 'NPV', 'ROC_AUC', 'PRC_AUC', 'MCC'
]

EXTENDED_METRICS = [
    'Kappa', 'Brier', 'ECE', 'LR_Plus', 'LR_Minus', 'DOR', 'Youden_J'
]

ALL_METRICS = CORE_METRICS + EXTENDED_METRICS

PLOT_METRICS = [
    'Accuracy', 'Balanced_Accuracy', 'F1', 'Precision', 'Recall',
    'Specificity', 'NPV', 'ROC_AUC', 'PRC_AUC'
]

# ---------------------------------------------------------------------------
# PATIENT-LEVEL TRACKING
# ---------------------------------------------------------------------------
TRACK_PATIENT_PREDICTIONS = True
STABILITY_THRESHOLD = 0.7

# ---------------------------------------------------------------------------
# ROC/PRC CURVE PARAMETERS
# ---------------------------------------------------------------------------
ROC_N_BOOTSTRAP_CI = 1000

# ---------------------------------------------------------------------------
# FEATURE ANALYSIS PARAMETERS
# ---------------------------------------------------------------------------
FEATURE_CORRELATION_METHOD = 'spearman'
IMPORTANCE_N_PERMUTATIONS = 100

# ---------------------------------------------------------------------------
# LEARNING CURVE PARAMETERS
# ---------------------------------------------------------------------------
LEARNING_CURVE_CV_FOLDS = 5



# =============================================================================
# IMPORTS
# =============================================================================

import sys
import json
import pickle
import warnings
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
import matplotlib.transforms as transforms
import seaborn as sns
from scipy import stats
from scipy.special import ndtri

from sklearn.model_selection import (
    RepeatedStratifiedKFold, StratifiedKFold,
    cross_val_predict, permutation_test_score
)
from sklearn.preprocessing import StandardScaler
from sklearn.naive_bayes import GaussianNB
from sklearn.dummy import DummyClassifier
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.pipeline import Pipeline
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score,
    precision_score, recall_score, roc_auc_score,
    average_precision_score, confusion_matrix,
    roc_curve, precision_recall_curve, matthews_corrcoef,
    cohen_kappa_score, brier_score_loss
)
from sklearn.base import clone, BaseEstimator, ClassifierMixin
from sklearn.utils import resample

warnings.filterwarnings('ignore')

try:
    import shap
    HAS_SHAP = True
except ImportError:
    HAS_SHAP = False
    print("[WARNING] SHAP not installed. SHAP analysis will be skipped.")

try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False
    print("[WARNING] lifelines not installed. KM analysis will be limited.")

try:
    from tqdm.auto import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    def tqdm(x, **kwargs): return x

# Publication-quality plot defaults with TrueType fonts for Adobe Illustrator
plt.rcParams.update({
    'font.size': FONT_SIZE,
    'font.family': 'sans-serif',
    'axes.labelsize': LABEL_SIZE,
    'axes.titlesize': TITLE_SIZE,
    'xtick.labelsize': TICK_SIZE,
    'ytick.labelsize': TICK_SIZE,
    'legend.fontsize': TICK_SIZE,
    'figure.titlesize': TITLE_SIZE,
    'figure.dpi': 100,
    'savefig.dpi': FIGURE_DPI,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none'
})

# Define numpy-dependent configuration variables (after imports)
ROC_MEAN_FPR = np.linspace(0, 1, 100)
PRC_MEAN_RECALL = np.linspace(0, 1, 100)
LEARNING_CURVE_TRAIN_SIZES = np.linspace(0.2, 1.0, 10)



# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

LOG_PATH = os.path.join(OUTPUT_DIR, "analysis_log.txt")


def log(msg: str, level: str = "INFO"):
    """Log message to console and file."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] [{level}] {msg}"
    if VERBOSE:
        try:
            print(line)
        except UnicodeEncodeError:
            print(line.encode('ascii', 'replace').decode('ascii'))
    try:
        with open(LOG_PATH, "a", encoding='utf-8') as f:
            f.write(line + "\n")
    except:
        pass


def savefig(name: str, fig=None, subdir: str = None):
    """Save figure in multiple formats."""
    if fig is None:
        fig = plt.gcf()

    save_dir = OUTPUT_DIR if subdir is None else os.path.join(OUTPUT_DIR, subdir)
    os.makedirs(save_dir, exist_ok=True)

    for fmt in SAVE_FORMATS:
        path = os.path.join(save_dir, f"{name}.{fmt}")
        fig.savefig(path, bbox_inches='tight', dpi=FIGURE_DPI, facecolor='white')

    plt.close(fig)
    log(f"Saved figure: {name}")


def save_excel(path: str, sheets: Dict[str, pd.DataFrame]):
    """Save multiple DataFrames to Excel."""
    try:
        with pd.ExcelWriter(path, engine='openpyxl') as writer:
            for name, df in sheets.items():
                if isinstance(df, dict):
                    df = pd.DataFrame(df)
                df.to_excel(writer, sheet_name=str(name)[:31], index=False)
        log(f"Saved Excel: {os.path.basename(path)}")
        return path
    except Exception as e:
        log(f"Excel save failed: {e}", "ERROR")
        return None


def create_configuration_summary():
    """Create a summary of all configuration parameters."""
    config = {
        'Parameter': [],
        'Value': [],
        'Description': []
    }

    config['Parameter'].extend(['N_BOOTSTRAP', 'CV_N_REPEATS', 'CV_N_FOLDS', 'TOTAL_EVALUATIONS'])
    config['Value'].extend([N_BOOTSTRAP, CV_N_REPEATS, CV_N_FOLDS, TOTAL_EVALUATIONS])
    config['Description'].extend(['Bootstrap iterations', 'CV repetitions', 'CV folds', 'Total evaluations'])

    config['Parameter'].extend(['CI_LEVEL', 'CI_METHOD', 'N_PERMUTATIONS'])
    config['Value'].extend([CI_LEVEL, 'percentile (direct from bootstrap)', N_PERMUTATIONS])
    config['Description'].extend(['CI level', 'CI method', 'Permutation tests'])

    config['Parameter'].extend(['NB_FEATURES', 'COX_FEATURES'])
    config['Value'].extend([str(NB_FEATURES), str(COX_FEATURES)])
    config['Description'].extend(['ML model features', 'Cox model features'])

    for feat, beta in COX_BETAS.items():
        config['Parameter'].append(f'COX_BETA_{feat}')
        config['Value'].append(beta)
        config['Description'].append(f'Cox beta coefficient for {feat}')

    config['Parameter'].extend(['RANDOM_STATE', 'STABILITY_THRESHOLD', 'RISK_THRESHOLD_METHOD'])
    config['Value'].extend([RANDOM_STATE, STABILITY_THRESHOLD, RISK_THRESHOLD_METHOD])
    config['Description'].extend(['Random seed', 'Stability threshold', 'Risk threshold method'])

    return pd.DataFrame(config)



# =============================================================================
# DELONG TEST FOR AUC COMPARISON
# =============================================================================

def compute_midrank(x):
    """Compute midranks for DeLong test."""
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=np.float64)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1)
        i = j
    T2 = np.empty(N, dtype=np.float64)
    T2[J] = T + 1
    return T2


def fastDeLong(predictions_sorted_transposed, label_1_count):
    """Fast DeLong AUC computation."""
    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m
    positive_examples = predictions_sorted_transposed[:, :m]
    negative_examples = predictions_sorted_transposed[:, m:]
    k = predictions_sorted_transposed.shape[0]

    aucs = np.zeros(k)
    for r in range(k):
        midrank = compute_midrank(predictions_sorted_transposed[r, :])
        aucs[r] = (np.sum(midrank[:m]) - m * (m + 1) / 2) / (m * n)

    return aucs


def delong_test(y_true: np.ndarray, y_pred1: np.ndarray, y_pred2: np.ndarray) -> Dict:
    """
    Perform DeLong test for comparing two ROC AUCs.

    Reference: DeLong et al. (1988)

    Returns dict with AUCs, difference, z-score, and p-value.
    """
    y_true = np.asarray(y_true)
    y_pred1 = np.asarray(y_pred1)
    y_pred2 = np.asarray(y_pred2)

    auc1 = roc_auc_score(y_true, y_pred1)
    auc2 = roc_auc_score(y_true, y_pred2)

    order = np.argsort(y_true)[::-1]
    y_true_sorted = y_true[order]
    y_pred1_sorted = y_pred1[order]
    y_pred2_sorted = y_pred2[order]

    label_1_count = int(np.sum(y_true))
    predictions_sorted = np.vstack([y_pred1_sorted, y_pred2_sorted])
    aucs = fastDeLong(predictions_sorted, label_1_count)

    # Compute covariance matrix using structural components
    m = label_1_count
    n = len(y_true) - m

    positive_examples = predictions_sorted[:, :m]
    negative_examples = predictions_sorted[:, m:]

    # V10 and V01 placement value matrices
    V10 = np.zeros((2, m))
    V01 = np.zeros((2, n))

    for i in range(2):
        for j in range(m):
            V10[i, j] = np.mean(positive_examples[i, j] > negative_examples[i, :]) + \
                       0.5 * np.mean(positive_examples[i, j] == negative_examples[i, :])
        for j in range(n):
            V01[i, j] = np.mean(negative_examples[i, j] < positive_examples[i, :]) + \
                       0.5 * np.mean(negative_examples[i, j] == positive_examples[i, :])

    S10 = np.cov(V10)
    S01 = np.cov(V01)
    S = S10 / m + S01 / n

    # DeLong test statistic
    L = np.array([1, -1])
    var_diff = L @ S @ L.T

    if var_diff <= 0:
        z_score = 0
        p_value = 1.0
    else:
        z_score = (auc1 - auc2) / np.sqrt(var_diff)
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

    return {
        'auc1': auc1,
        'auc2': auc2,
        'difference': auc1 - auc2,
        'z_score': z_score,
        'p_value': p_value,
        'significant': p_value < DELONG_ALPHA
    }



# =============================================================================
# NET RECLASSIFICATION IMPROVEMENT (NRI) & IDI
# =============================================================================

def compute_nri_idi(y_true: np.ndarray, y_proba_old: np.ndarray,
                    y_proba_new: np.ndarray, threshold: float = 0.5) -> Dict:
    """
    Compute Net Reclassification Improvement (NRI) and
    Integrated Discrimination Improvement (IDI).

    NRI measures how well a new model reclassifies subjects.
    IDI measures improvement in discrimination slopes.

    Reference: Pencina et al. (2008)

    Parameters:
    -----------
    y_true : array-like
        True binary labels
    y_proba_old : array-like
        Predicted probabilities from old/reference model
    y_proba_new : array-like
        Predicted probabilities from new model
    threshold : float
        Threshold for categorical NRI

    Returns:
    --------
    dict with NRI, IDI, and components
    """
    y_true = np.asarray(y_true)
    y_proba_old = np.asarray(y_proba_old)
    y_proba_new = np.asarray(y_proba_new)

    n = len(y_true)
    events = y_true == 1
    non_events = y_true == 0
    n_events = np.sum(events)
    n_non_events = np.sum(non_events)

    pred_old = (y_proba_old >= threshold).astype(int)
    pred_new = (y_proba_new >= threshold).astype(int)

    # Categorical NRI: for events, reclassification up is good
    events_up = np.sum((pred_new > pred_old) & events)
    events_down = np.sum((pred_new < pred_old) & events)

    # For non-events, reclassification down is good
    non_events_up = np.sum((pred_new > pred_old) & non_events)
    non_events_down = np.sum((pred_new < pred_old) & non_events)

    nri_events = (events_up - events_down) / n_events if n_events > 0 else 0
    nri_non_events = (non_events_down - non_events_up) / n_non_events if n_non_events > 0 else 0
    nri = nri_events + nri_non_events

    se_nri_events = np.sqrt((events_up + events_down) / n_events**2) if n_events > 0 else 0
    se_nri_non_events = np.sqrt((non_events_up + non_events_down) / n_non_events**2) if n_non_events > 0 else 0
    se_nri = np.sqrt(se_nri_events**2 + se_nri_non_events**2)

    z_nri = nri / se_nri if se_nri > 0 else 0
    p_nri = 2 * (1 - stats.norm.cdf(abs(z_nri)))

    # Continuous NRI
    prob_change_events = y_proba_new[events] - y_proba_old[events]
    prob_change_non_events = y_proba_new[non_events] - y_proba_old[non_events]

    nri_cont_events = np.mean(prob_change_events > 0) - np.mean(prob_change_events < 0) if n_events > 0 else 0
    nri_cont_non_events = np.mean(prob_change_non_events < 0) - np.mean(prob_change_non_events > 0) if n_non_events > 0 else 0
    nri_continuous = nri_cont_events + nri_cont_non_events

    # IDI: discrimination slopes
    slope_old = np.mean(y_proba_old[events]) - np.mean(y_proba_old[non_events]) if n_events > 0 and n_non_events > 0 else 0
    slope_new = np.mean(y_proba_new[events]) - np.mean(y_proba_new[non_events]) if n_events > 0 and n_non_events > 0 else 0
    idi = slope_new - slope_old

    var_new_events = np.var(y_proba_new[events]) if n_events > 1 else 0
    var_new_non = np.var(y_proba_new[non_events]) if n_non_events > 1 else 0
    var_old_events = np.var(y_proba_old[events]) if n_events > 1 else 0
    var_old_non = np.var(y_proba_old[non_events]) if n_non_events > 1 else 0

    se_idi = np.sqrt(var_new_events/n_events + var_new_non/n_non_events +
                     var_old_events/n_events + var_old_non/n_non_events)

    z_idi = idi / se_idi if se_idi > 0 else 0
    p_idi = 2 * (1 - stats.norm.cdf(abs(z_idi)))

    return {
        'NRI': nri, 'NRI_events': nri_events, 'NRI_non_events': nri_non_events,
        'NRI_SE': se_nri, 'NRI_z': z_nri, 'NRI_p': p_nri,
        'NRI_continuous': nri_continuous, 'NRI_cont_events': nri_cont_events,
        'NRI_cont_non_events': nri_cont_non_events,
        'IDI': idi, 'IDI_SE': se_idi, 'IDI_z': z_idi, 'IDI_p': p_idi,
        'slope_old': slope_old, 'slope_new': slope_new,
        'events_up': events_up, 'events_down': events_down,
        'non_events_up': non_events_up, 'non_events_down': non_events_down
    }



# =============================================================================
# DECISION CURVE ANALYSIS (DCA)
# =============================================================================

def decision_curve_analysis(y_true: np.ndarray, y_proba: np.ndarray,
                            thresholds: np.ndarray = None) -> pd.DataFrame:
    """
    Perform Decision Curve Analysis.

    Calculates net benefit at different threshold probabilities.

    Reference: Vickers & Elkin (2006)

    Net Benefit = (TP/n) - (FP/n) * (pt / (1-pt))
    where pt is threshold probability
    """
    y_true = np.asarray(y_true)
    y_proba = np.asarray(y_proba)
    n = len(y_true)

    if thresholds is None:
        thresholds = np.linspace(DCA_THRESHOLD_RANGE[0], DCA_THRESHOLD_RANGE[1], DCA_N_THRESHOLDS)

    results = []
    for pt in thresholds:
        y_pred = (y_proba >= pt).astype(int)
        tp = np.sum((y_pred == 1) & (y_true == 1))
        fp = np.sum((y_pred == 1) & (y_true == 0))

        if pt < 1:
            net_benefit = (tp / n) - (fp / n) * (pt / (1 - pt))
        else:
            net_benefit = 0

        prevalence = np.mean(y_true)
        treat_all_nb = prevalence - (1 - prevalence) * (pt / (1 - pt)) if pt < 1 else 0

        results.append({
            'threshold': pt, 'net_benefit': net_benefit,
            'treat_all': treat_all_nb, 'treat_none': 0, 'tp': tp, 'fp': fp
        })

    return pd.DataFrame(results)


# =============================================================================
# BRIER SCORE DECOMPOSITION
# =============================================================================

def brier_score_decomposition(y_true: np.ndarray, y_proba: np.ndarray,
                               n_bins: int = 10) -> Dict:
    """
    Decompose Brier score into reliability, resolution, and uncertainty.

    Brier = Reliability - Resolution + Uncertainty

    - Reliability (calibration): How well probabilities match frequencies
    - Resolution: How much predictions vary from base rate
    - Uncertainty: Inherent uncertainty in the outcome
    """
    y_true = np.asarray(y_true)
    y_proba = np.asarray(y_proba)
    n = len(y_true)

    brier = brier_score_loss(y_true, y_proba)
    base_rate = np.mean(y_true)
    uncertainty = base_rate * (1 - base_rate)

    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_indices = np.digitize(y_proba, bin_edges[1:-1])

    reliability = 0
    resolution = 0

    for k in range(n_bins):
        mask = bin_indices == k
        n_k = np.sum(mask)
        if n_k > 0:
            f_k = np.mean(y_proba[mask])
            o_k = np.mean(y_true[mask])
            reliability += n_k * (f_k - o_k) ** 2
            resolution += n_k * (o_k - base_rate) ** 2

    reliability /= n
    resolution /= n

    # Brier skill score (relative to climatology)
    brier_skill = 1 - brier / uncertainty if uncertainty > 0 else 0

    return {
        'brier': brier, 'reliability': reliability, 'resolution': resolution,
        'uncertainty': uncertainty, 'brier_skill_score': brier_skill,
        'calibration_component': reliability, 'discrimination_component': resolution
    }


# =============================================================================
# THRESHOLD OPTIMIZATION
# =============================================================================

def optimize_threshold(y_true: np.ndarray, y_proba: np.ndarray,
                       method: str = 'youden') -> Dict:
    """
    Find optimal classification threshold.

    Methods:
    - 'youden': Maximize Youden's J (sensitivity + specificity - 1)
    - 'f1': Maximize F1 score
    - 'balanced': Minimize |sensitivity - specificity|
    """
    fpr, tpr, thresholds = roc_curve(y_true, y_proba)
    results = {}

    # Youden's J
    j_scores = tpr - fpr
    optimal_idx_youden = np.argmax(j_scores)
    results['youden'] = {
        'threshold': thresholds[optimal_idx_youden],
        'sensitivity': tpr[optimal_idx_youden],
        'specificity': 1 - fpr[optimal_idx_youden],
        'youden_j': j_scores[optimal_idx_youden]
    }

    # F1 optimization
    f1_scores = []
    for thresh in thresholds:
        y_pred = (y_proba >= thresh).astype(int)
        f1_scores.append(f1_score(y_true, y_pred, zero_division=0))
    optimal_idx_f1 = np.argmax(f1_scores)
    results['f1'] = {
        'threshold': thresholds[optimal_idx_f1],
        'f1_score': f1_scores[optimal_idx_f1]
    }

    # Balanced (equal sensitivity/specificity)
    sens_spec_diff = np.abs(tpr - (1 - fpr))
    optimal_idx_balanced = np.argmin(sens_spec_diff)
    results['balanced'] = {
        'threshold': thresholds[optimal_idx_balanced],
        'sensitivity': tpr[optimal_idx_balanced],
        'specificity': 1 - fpr[optimal_idx_balanced]
    }

    return results



# =============================================================================
# COX MODEL AS SKLEARN-COMPATIBLE CLASSIFIER
# =============================================================================

class CoxRiskClassifier(BaseEstimator, ClassifierMixin):
    """
    Wrapper to use Cox proportional hazards model as a binary classifier.
    Uses pre-defined beta coefficients to compute linear predictor (risk score).
    """

    def __init__(self, betas: Dict[str, float], features: List[str],
                 threshold_method: str = 'median', fixed_threshold: float = 0.5):
        self.betas = betas
        self.features = features
        self.threshold_method = threshold_method
        self.fixed_threshold = fixed_threshold
        self.threshold_ = None
        self.scaler_ = None
        self.classes_ = np.array([0, 1])

    def _compute_risk_score(self, X: np.ndarray) -> np.ndarray:
        """Compute linear predictor (sum of beta * x)."""
        risk = np.zeros(X.shape[0])
        for i, feat in enumerate(self.features):
            if feat in self.betas:
                risk += self.betas[feat] * X[:, i]
        return risk

    def _risk_to_probability(self, risk: np.ndarray) -> np.ndarray:
        """Convert risk score to probability using sigmoid transformation."""
        centered_risk = risk - np.mean(risk) if len(risk) > 1 else risk
        prob = 1 / (1 + np.exp(-centered_risk))
        return prob

    def fit(self, X: np.ndarray, y: np.ndarray):
        """Fit the classifier (determine threshold from training data)."""
        self.scaler_ = StandardScaler()
        X_scaled = self.scaler_.fit_transform(X)
        risk_scores = self._compute_risk_score(X_scaled)

        if self.threshold_method == 'median':
            self.threshold_ = np.median(risk_scores)
        elif self.threshold_method == 'youden':
            probs = self._risk_to_probability(risk_scores)
            fpr, tpr, thresholds_roc = roc_curve(y, probs)
            j_scores = tpr - fpr
            optimal_idx = np.argmax(j_scores)
            self.threshold_ = np.percentile(risk_scores,
                                           (1 - thresholds_roc[optimal_idx]) * 100)
        else:
            self.threshold_ = np.percentile(risk_scores, self.fixed_threshold * 100)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict class labels."""
        X_scaled = self.scaler_.transform(X)
        risk_scores = self._compute_risk_score(X_scaled)
        return (risk_scores > self.threshold_).astype(int)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """Predict class probabilities."""
        X_scaled = self.scaler_.transform(X)
        risk_scores = self._compute_risk_score(X_scaled)
        probs = self._risk_to_probability(risk_scores)
        return np.column_stack([1 - probs, probs])

    def get_risk_scores(self, X: np.ndarray) -> np.ndarray:
        """Get raw risk scores (linear predictor)."""
        X_scaled = self.scaler_.transform(X)
        return self._compute_risk_score(X_scaled)


# =============================================================================
# NAIVE BAYES CLASSIFIER WRAPPER
# =============================================================================

class NBChampionClassifier(BaseEstimator, ClassifierMixin):
    """Wrapper for the champion NB model that handles scaling internally."""

    def __init__(self, features: List[str], calibrate: bool = True,
                 calibration_cv: int = 3):
        self.features = features
        self.calibrate = calibrate
        self.calibration_cv = calibration_cv
        self.model_ = None
        self.scaler_ = None
        self.classes_ = np.array([0, 1])

    def fit(self, X: np.ndarray, y: np.ndarray):
        """Fit the NB classifier."""
        self.scaler_ = StandardScaler()
        X_scaled = self.scaler_.fit_transform(X)
        base_clf = GaussianNB()

        if self.calibrate and len(np.unique(y)) > 1:
            min_class_count = min(np.sum(y == 0), np.sum(y == 1))
            cv_folds = min(self.calibration_cv, min_class_count)
            if cv_folds >= 2:
                self.model_ = CalibratedClassifierCV(
                    base_clf, method='sigmoid', cv=cv_folds
                )
            else:
                self.model_ = base_clf
        else:
            self.model_ = base_clf

        self.model_.fit(X_scaled, y)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict class labels."""
        X_scaled = self.scaler_.transform(X)
        return self.model_.predict(X_scaled)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """Predict class probabilities."""
        X_scaled = self.scaler_.transform(X)
        return self.model_.predict_proba(X_scaled)



# =============================================================================
# METRICS COMPUTATION
# =============================================================================

def compute_ece(y_true: np.ndarray, y_proba: np.ndarray, n_bins: int = 10) -> float:
    """Compute Expected Calibration Error."""
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    total_samples = len(y_true)

    for i in range(n_bins):
        bin_lower, bin_upper = bin_boundaries[i], bin_boundaries[i + 1]
        in_bin = (y_proba >= bin_lower) & (y_proba < bin_upper)
        if np.sum(in_bin) > 0:
            bin_accuracy = np.mean(y_true[in_bin])
            bin_confidence = np.mean(y_proba[in_bin])
            bin_size = np.sum(in_bin)
            ece += (bin_size / total_samples) * np.abs(bin_accuracy - bin_confidence)

    return ece


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray,
                   y_proba: np.ndarray = None) -> Dict[str, float]:
    """Compute comprehensive classification metrics."""
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    metrics = {
        'Accuracy': accuracy_score(y_true, y_pred),
        'Balanced_Accuracy': balanced_accuracy_score(y_true, y_pred),
        'F1': f1_score(y_true, y_pred, zero_division=0),
        'Precision': precision_score(y_true, y_pred, zero_division=0),
        'Recall': recall_score(y_true, y_pred, zero_division=0),
        'Specificity': tn / (tn + fp) if (tn + fp) > 0 else 0,
        'NPV': tn / (tn + fn) if (tn + fn) > 0 else 0,
        'MCC': matthews_corrcoef(y_true, y_pred),
        'Kappa': cohen_kappa_score(y_true, y_pred),
        'TP': tp, 'TN': tn, 'FP': fp, 'FN': fn
    }

    if y_proba is not None:
        try:
            metrics['ROC_AUC'] = roc_auc_score(y_true, y_proba)
        except:
            metrics['ROC_AUC'] = 0.5
        try:
            metrics['PRC_AUC'] = average_precision_score(y_true, y_proba)
        except:
            metrics['PRC_AUC'] = np.mean(y_true)

        metrics['Brier'] = brier_score_loss(y_true, y_proba)
        metrics['ECE'] = compute_ece(y_true, y_proba)

        sens = metrics['Recall']
        spec = metrics['Specificity']
        metrics['LR_Plus'] = sens / (1 - spec) if (1 - spec) > 0.001 else np.inf
        metrics['LR_Minus'] = (1 - sens) / spec if spec > 0.001 else np.inf

    return metrics


def bca_ci(data: np.ndarray, stat_func, n_bootstrap: int = 5000,
           ci_level: float = 0.95, seed: int = 42) -> Tuple[float, float, float]:
    """Compute BCa bootstrap confidence interval."""
    rng = np.random.RandomState(seed)
    data = np.asarray(data)
    n = len(data)

    if n == 0:
        return (np.nan, np.nan, np.nan)
    try:
        theta_hat = stat_func(data)
    except:
        return (np.nan, np.nan, np.nan)

    theta_boot = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        boot_sample = data[rng.randint(0, n, size=n)]
        try:
            theta_boot[i] = stat_func(boot_sample)
        except:
            theta_boot[i] = np.nan

    theta_boot = theta_boot[~np.isnan(theta_boot)]
    if len(theta_boot) < 100:
        return (theta_hat, np.nan, np.nan)

    prop_less = np.mean(theta_boot < theta_hat)
    if prop_less == 0:
        prop_less = 1 / (2 * n_bootstrap)
    elif prop_less == 1:
        prop_less = 1 - 1 / (2 * n_bootstrap)
    z0 = ndtri(prop_less)

    theta_jack = np.zeros(n)
    for i in range(n):
        jack_sample = np.delete(data, i)
        try:
            theta_jack[i] = stat_func(jack_sample)
        except:
            theta_jack[i] = np.nan

    theta_jack = theta_jack[~np.isnan(theta_jack)]
    if len(theta_jack) < 2:
        alpha = 1 - ci_level
        ci_lower = np.percentile(theta_boot, alpha/2 * 100)
        ci_upper = np.percentile(theta_boot, (1 - alpha/2) * 100)
        return (theta_hat, ci_lower, ci_upper)

    theta_jack_mean = np.mean(theta_jack)
    num = np.sum((theta_jack_mean - theta_jack) ** 3)
    denom = 6 * (np.sum((theta_jack_mean - theta_jack) ** 2) ** 1.5)

    if denom == 0:
        a = 0
    else:
        a = num / denom

    alpha = 1 - ci_level
    z_alpha_low = ndtri(alpha / 2)
    z_alpha_high = ndtri(1 - alpha / 2)

    def adjusted_percentile(z_alpha):
        numer = z0 + z_alpha
        denom_adj = 1 - a * numer
        if denom_adj == 0:
            return 0.5
        adjusted_z = z0 + numer / denom_adj
        return stats.norm.cdf(adjusted_z)

    p_low = adjusted_percentile(z_alpha_low) * 100
    p_high = adjusted_percentile(z_alpha_high) * 100
    p_low = np.clip(p_low, 0.1, 99.9)
    p_high = np.clip(p_high, 0.1, 99.9)

    if np.isnan(p_low) or np.isnan(p_high) or p_low >= p_high:
        alpha = 1 - ci_level
        ci_lower = np.percentile(theta_boot, alpha/2 * 100)
        ci_upper = np.percentile(theta_boot, (1 - alpha/2) * 100)
        return (theta_hat, ci_lower, ci_upper)

    ci_lower = np.percentile(theta_boot, p_low)
    ci_upper = np.percentile(theta_boot, p_high)
    return (theta_hat, ci_lower, ci_upper)


def percentile_ci(data: np.ndarray, ci_level: float = 0.95) -> Tuple[float, float, float]:
    """
    Compute percentile confidence interval directly from bootstrap distribution.

    When you already have bootstrap samples (e.g., 200 performance metrics from
    bootstrap CV), the distribution IS the bootstrap distribution. No secondary
    bootstrap is needed; just take percentiles directly.

    Args:
        data: Array of bootstrap samples (e.g., 200 F1 scores from bootstrap CV)
        ci_level: Confidence level (default 0.95 for 95% CI)

    Returns:
        Tuple of (mean, ci_lower, ci_upper)
    """
    data = np.asarray(data)
    data = data[~np.isnan(data)]

    if len(data) == 0:
        return (np.nan, np.nan, np.nan)

    alpha = 1 - ci_level
    mean_val = np.mean(data)
    ci_lower = np.percentile(data, alpha / 2 * 100)
    ci_upper = np.percentile(data, (1 - alpha / 2) * 100)

    return (mean_val, ci_lower, ci_upper)



# =============================================================================
# BOOTSTRAP CV EVALUATION
# =============================================================================

def run_bootstrap_cv_evaluation(X_nb: np.ndarray, X_cox: np.ndarray, y: np.ndarray,
                                 patient_ids: np.ndarray = None,
                                 n_bootstrap: int = 200,
                                 n_repeats: int = 20,
                                 n_folds: int = 5,
                                 seed: int = 42) -> Dict:
    """Run bootstrap CV evaluation for all three models: NB, Cox, ZeroR.

    Includes patient-level prediction tracking for stability analysis.
    """
    log(f"\nRunning Bootstrap CV: {n_bootstrap} bootstrap x {n_repeats}x{n_folds} CV")
    log(f"Total evaluations: {n_bootstrap * n_repeats * n_folds}")

    rng = np.random.RandomState(seed)
    n_samples = len(y)

    if patient_ids is None:
        patient_ids = np.arange(n_samples)

    all_metrics_nb = defaultdict(list)
    all_metrics_cox = defaultdict(list)
    all_metrics_zeror = defaultdict(list)

    all_y_true = []
    all_y_pred_nb = []
    all_y_pred_cox = []
    all_y_pred_zeror = []
    all_y_proba_nb = []
    all_y_proba_cox = []

    patient_predictions_nb = defaultdict(list)
    patient_predictions_cox = defaultdict(list)
    patient_probas_nb = defaultdict(list)
    patient_probas_cox = defaultdict(list)
    patient_true = {}

    bootstrap_roc_nb = []
    bootstrap_roc_cox = []
    bootstrap_prc_nb = []
    bootstrap_prc_cox = []

    for boot_idx in tqdm(range(n_bootstrap), desc="Bootstrap iterations"):
        boot_indices = rng.choice(n_samples, size=n_samples, replace=True)

        X_nb_boot = X_nb[boot_indices]
        X_cox_boot = X_cox[boot_indices]
        y_boot = y[boot_indices]
        patient_ids_boot = patient_ids[boot_indices]

        cv = RepeatedStratifiedKFold(n_splits=n_folds, n_repeats=n_repeats,
                                     random_state=seed + boot_idx)

        fold_metrics_nb = defaultdict(list)
        fold_metrics_cox = defaultdict(list)
        fold_metrics_zeror = defaultdict(list)

        boot_y_true = []
        boot_y_proba_nb = []
        boot_y_proba_cox = []

        for train_idx, val_idx in cv.split(X_nb_boot, y_boot):
            X_nb_tr, X_nb_val = X_nb_boot[train_idx], X_nb_boot[val_idx]
            X_cox_tr, X_cox_val = X_cox_boot[train_idx], X_cox_boot[val_idx]
            y_tr, y_val = y_boot[train_idx], y_boot[val_idx]
            val_patient_ids = patient_ids_boot[val_idx]

            if len(np.unique(y_val)) < 2:
                continue

            try:
                nb_clf = NBChampionClassifier(features=NB_FEATURES, calibrate=True)
                nb_clf.fit(X_nb_tr, y_tr)
                y_pred_nb = nb_clf.predict(X_nb_val)
                y_proba_nb = nb_clf.predict_proba(X_nb_val)[:, 1]
            except:
                y_pred_nb = np.zeros(len(y_val))
                y_proba_nb = np.full(len(y_val), 0.5)

            try:
                cox_clf = CoxRiskClassifier(betas=COX_BETAS, features=COX_FEATURES,
                                           threshold_method=RISK_THRESHOLD_METHOD)
                cox_clf.fit(X_cox_tr, y_tr)
                y_pred_cox = cox_clf.predict(X_cox_val)
                y_proba_cox = cox_clf.predict_proba(X_cox_val)[:, 1]
            except:
                y_pred_cox = np.zeros(len(y_val))
                y_proba_cox = np.full(len(y_val), 0.5)

            zeror = DummyClassifier(strategy="most_frequent")
            zeror.fit(X_nb_tr, y_tr)
            y_pred_zeror = zeror.predict(X_nb_val)
            majority_prob = np.mean(y_tr)
            y_proba_zeror = np.full(len(y_val), majority_prob)

            metrics_nb = compute_metrics(y_val, y_pred_nb, y_proba_nb)
            metrics_cox = compute_metrics(y_val, y_pred_cox, y_proba_cox)
            metrics_zeror = compute_metrics(y_val, y_pred_zeror, y_proba_zeror)

            for metric, val in metrics_nb.items():
                fold_metrics_nb[metric].append(val)
            for metric, val in metrics_cox.items():
                fold_metrics_cox[metric].append(val)
            for metric, val in metrics_zeror.items():
                fold_metrics_zeror[metric].append(val)

            all_y_true.extend(y_val)
            all_y_pred_nb.extend(y_pred_nb)
            all_y_pred_cox.extend(y_pred_cox)
            all_y_pred_zeror.extend(y_pred_zeror)
            all_y_proba_nb.extend(y_proba_nb)
            all_y_proba_cox.extend(y_proba_cox)

            if TRACK_PATIENT_PREDICTIONS:
                for i, pid in enumerate(val_patient_ids):
                    patient_predictions_nb[pid].append(y_pred_nb[i])
                    patient_predictions_cox[pid].append(y_pred_cox[i])
                    patient_probas_nb[pid].append(y_proba_nb[i])
                    patient_probas_cox[pid].append(y_proba_cox[i])
                    patient_true[pid] = y_val[i]

            boot_y_true.extend(y_val)
            boot_y_proba_nb.extend(y_proba_nb)
            boot_y_proba_cox.extend(y_proba_cox)

        for metric in fold_metrics_nb:
            if fold_metrics_nb[metric]:
                all_metrics_nb[metric].append(np.mean(fold_metrics_nb[metric]))
        for metric in fold_metrics_cox:
            if fold_metrics_cox[metric]:
                all_metrics_cox[metric].append(np.mean(fold_metrics_cox[metric]))
        for metric in fold_metrics_zeror:
            if fold_metrics_zeror[metric]:
                all_metrics_zeror[metric].append(np.mean(fold_metrics_zeror[metric]))

        if len(boot_y_true) > 10 and len(np.unique(boot_y_true)) == 2:
            try:
                fpr_nb, tpr_nb, _ = roc_curve(boot_y_true, boot_y_proba_nb)
                fpr_cox, tpr_cox, _ = roc_curve(boot_y_true, boot_y_proba_cox)
                bootstrap_roc_nb.append((fpr_nb, tpr_nb))
                bootstrap_roc_cox.append((fpr_cox, tpr_cox))
                prec_nb, rec_nb, _ = precision_recall_curve(boot_y_true, boot_y_proba_nb)
                prec_cox, rec_cox, _ = precision_recall_curve(boot_y_true, boot_y_proba_cox)
                bootstrap_prc_nb.append((rec_nb, prec_nb))
                bootstrap_prc_cox.append((rec_cox, prec_cox))
            except:
                pass

    df_nb = pd.DataFrame(dict(all_metrics_nb))
    df_cox = pd.DataFrame(dict(all_metrics_cox))
    df_zeror = pd.DataFrame(dict(all_metrics_zeror))

    summary_nb = {}
    summary_cox = {}
    summary_zeror = {}

    # Compute CIs directly from bootstrap distribution (no secondary bootstrap needed)
    # The 200 bootstrap samples already represent the sampling distribution
    for metric in df_nb.columns:
        if metric not in ['TP', 'TN', 'FP', 'FN']:
            summary_nb[metric] = percentile_ci(df_nb[metric].values, CI_LEVEL)
            summary_cox[metric] = percentile_ci(df_cox[metric].values, CI_LEVEL)
            summary_zeror[metric] = percentile_ci(df_zeror[metric].values, CI_LEVEL)

    paired_tests_nb_vs_zeror = {}
    paired_tests_cox_vs_zeror = {}
    paired_tests_nb_vs_cox = {}
    cohens_d_nb_vs_zeror = {}
    cohens_d_cox_vs_zeror = {}
    cohens_d_nb_vs_cox = {}

    for metric in CORE_METRICS + ['MCC']:
        if metric in df_nb.columns:
            t_stat, p_val = stats.ttest_rel(df_nb[metric], df_zeror[metric])
            paired_tests_nb_vs_zeror[metric] = {'t_stat': t_stat, 'p_value': p_val}
            diff = df_nb[metric] - df_zeror[metric]
            cohens_d_nb_vs_zeror[metric] = np.mean(diff) / np.std(diff) if np.std(diff) > 0 else 0

            t_stat, p_val = stats.ttest_rel(df_cox[metric], df_zeror[metric])
            paired_tests_cox_vs_zeror[metric] = {'t_stat': t_stat, 'p_value': p_val}
            diff = df_cox[metric] - df_zeror[metric]
            cohens_d_cox_vs_zeror[metric] = np.mean(diff) / np.std(diff) if np.std(diff) > 0 else 0

            t_stat, p_val = stats.ttest_rel(df_nb[metric], df_cox[metric])
            paired_tests_nb_vs_cox[metric] = {'t_stat': t_stat, 'p_value': p_val}
            diff = df_nb[metric] - df_cox[metric]
            cohens_d_nb_vs_cox[metric] = np.mean(diff) / np.std(diff) if np.std(diff) > 0 else 0

    # Compute patient stability metrics
    patient_stability = {}
    if TRACK_PATIENT_PREDICTIONS and patient_predictions_nb:
        for pid in patient_true.keys():
            preds_nb = patient_predictions_nb.get(pid, [])
            preds_cox = patient_predictions_cox.get(pid, [])

            if len(preds_nb) > 0:
                # Stability = proportion of predictions matching majority vote
                majority_nb = 1 if np.mean(preds_nb) >= 0.5 else 0
                majority_cox = 1 if np.mean(preds_cox) >= 0.5 else 0

                stability_nb = np.mean([p == majority_nb for p in preds_nb])
                stability_cox = np.mean([p == majority_cox for p in preds_cox])

                patient_stability[pid] = {
                    'true_label': patient_true[pid],
                    'n_predictions': len(preds_nb),
                    'mean_proba_nb': np.mean(patient_probas_nb.get(pid, [0.5])),
                    'std_proba_nb': np.std(patient_probas_nb.get(pid, [0.5])),
                    'mean_proba_cox': np.mean(patient_probas_cox.get(pid, [0.5])),
                    'std_proba_cox': np.std(patient_probas_cox.get(pid, [0.5])),
                    'majority_vote_nb': majority_nb,
                    'majority_vote_cox': majority_cox,
                    'stability_nb': stability_nb,
                    'stability_cox': stability_cox,
                    'is_stable_nb': stability_nb >= STABILITY_THRESHOLD,
                    'is_stable_cox': stability_cox >= STABILITY_THRESHOLD,
                    'correct_nb': majority_nb == patient_true[pid],
                    'correct_cox': majority_cox == patient_true[pid]
                }

    # Debug logging of final computed results

    return {
        'bootstrap_metrics_nb': df_nb,
        'bootstrap_metrics_cox': df_cox,
        'bootstrap_metrics_zeror': df_zeror,
        'summary_nb': summary_nb,
        'summary_cox': summary_cox,
        'summary_zeror': summary_zeror,
        'paired_tests_nb_vs_zeror': paired_tests_nb_vs_zeror,
        'paired_tests_cox_vs_zeror': paired_tests_cox_vs_zeror,
        'paired_tests_nb_vs_cox': paired_tests_nb_vs_cox,
        'cohens_d_nb_vs_zeror': cohens_d_nb_vs_zeror,
        'cohens_d_cox_vs_zeror': cohens_d_cox_vs_zeror,
        'cohens_d_nb_vs_cox': cohens_d_nb_vs_cox,
        'all_predictions': {
            'y_true': np.array(all_y_true),
            'y_pred_nb': np.array(all_y_pred_nb),
            'y_pred_cox': np.array(all_y_pred_cox),
            'y_proba_nb': np.array(all_y_proba_nb),
            'y_proba_cox': np.array(all_y_proba_cox)
        },
        # Patient-level tracking data
        'patient_tracking': {
            'patient_predictions_nb': dict(patient_predictions_nb),
            'patient_predictions_cox': dict(patient_predictions_cox),
            'patient_probas_nb': dict(patient_probas_nb),
            'patient_probas_cox': dict(patient_probas_cox),
            'patient_true': patient_true,
            'patient_stability': patient_stability
        },
        # ROC/PRC curve data for CI visualization
        'roc_curves': {
            'bootstrap_roc_nb': bootstrap_roc_nb,
            'bootstrap_roc_cox': bootstrap_roc_cox
        },
        'prc_curves': {
            'bootstrap_prc_nb': bootstrap_prc_nb,
            'bootstrap_prc_cox': bootstrap_prc_cox
        }
    }

# =============================================================================
# PERMUTATION TESTING
# =============================================================================

def run_permutation_test(X: np.ndarray, y: np.ndarray, model_class,
                         model_kwargs: Dict, n_permutations: int = 1000,
                         cv_folds: int = 5, seed: int = 42) -> Dict:
    """Run permutation test to assess statistical significance."""

    log(f"\nRunning permutation test ({n_permutations} permutations)...")

    rng = np.random.RandomState(seed)
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=seed)

    # Actual performance
    actual_metrics = defaultdict(list)

    for train_idx, val_idx in cv.split(X, y):
        X_tr, X_val = X[train_idx], X[val_idx]
        y_tr, y_val = y[train_idx], y[val_idx]

        model = model_class(**model_kwargs)
        model.fit(X_tr, y_tr)

        y_pred = model.predict(X_val)
        y_proba = model.predict_proba(X_val)[:, 1]

        metrics = compute_metrics(y_val, y_pred, y_proba)
        for m, v in metrics.items():
            actual_metrics[m].append(v)

    actual_means = {m: np.mean(v) for m, v in actual_metrics.items()}

    # Permutation distribution
    perm_distributions = {m: [] for m in actual_means}

    for i in tqdm(range(n_permutations), desc="Permutation test"):
        y_perm = rng.permutation(y)

        perm_metrics = defaultdict(list)
        for train_idx, val_idx in cv.split(X, y_perm):
            X_tr, X_val = X[train_idx], X[val_idx]
            y_tr, y_val = y_perm[train_idx], y_perm[val_idx]

            try:
                model = model_class(**model_kwargs)
                model.fit(X_tr, y_tr)
                y_pred = model.predict(X_val)
                y_proba = model.predict_proba(X_val)[:, 1]
                metrics = compute_metrics(y_val, y_pred, y_proba)
                for m, v in metrics.items():
                    perm_metrics[m].append(v)
            except:
                pass

        for m, v in perm_metrics.items():
            if v:
                perm_distributions[m].append(np.mean(v))

    # Compute p-values
    p_values = {}
    for m, actual in actual_means.items():
        if m in perm_distributions and perm_distributions[m]:
            perm_vals = np.array(perm_distributions[m])
            p_values[m] = np.mean(perm_vals >= actual)

    return {
        'actual_metrics': actual_means,
        'p_values': p_values,
        'perm_distributions': perm_distributions
    }


# =============================================================================
# RISK STRATIFICATION ANALYSIS
# =============================================================================

def run_risk_stratification(X_nb: np.ndarray, X_cox: np.ndarray, y: np.ndarray,
                            time: np.ndarray, event: np.ndarray,
                            patient_ids: np.ndarray) -> Dict:
    """Perform risk stratification and Kaplan-Meier analysis."""

    log("\n" + "="*70)
    log("RISK STRATIFICATION ANALYSIS")
    log("="*70)

    nb_clf = NBChampionClassifier(features=NB_FEATURES, calibrate=True)
    nb_clf.fit(X_nb, y)

    cox_clf = CoxRiskClassifier(betas=COX_BETAS, features=COX_FEATURES,
                                threshold_method=RISK_THRESHOLD_METHOD)
    cox_clf.fit(X_cox, y)

    y_pred_nb = nb_clf.predict(X_nb)
    y_proba_nb = nb_clf.predict_proba(X_nb)[:, 1]

    y_pred_cox = cox_clf.predict(X_cox)
    y_proba_cox = cox_clf.predict_proba(X_cox)[:, 1]
    risk_scores_cox = cox_clf.get_risk_scores(X_cox)

    risk_group_nb = np.where(y_proba_nb >= 0.5, 'High', 'Low')
    risk_threshold_cox = np.median(risk_scores_cox)
    risk_group_cox = np.where(risk_scores_cox >= risk_threshold_cox, 'High', 'Low')

    df_risk = pd.DataFrame({
        'Patient_ID': patient_ids,
        'True_Label': ['S' if yi == 1 else 'L' for yi in y],
        'PFS_days': time,
        'Event': event,
        'NB_Probability': y_proba_nb,
        'NB_Prediction': ['S' if yi == 1 else 'L' for yi in y_pred_nb],
        'NB_Risk_Group': risk_group_nb,
        'Cox_Risk_Score': risk_scores_cox,
        'Cox_Probability': y_proba_cox,
        'Cox_Prediction': ['S' if yi == 1 else 'L' for yi in y_pred_cox],
        'Cox_Risk_Group': risk_group_cox
    })

    agreement = np.mean(y_pred_nb == y_pred_cox)
    kappa = cohen_kappa_score(y_pred_nb, y_pred_cox)

    log(f"Model agreement: {agreement:.1%}")
    log(f"Cohen's Kappa (NB vs Cox): {kappa:.3f}")

    # Agreement characterization
    agreement_analysis = characterize_disagreements(y, y_pred_nb, y_pred_cox,
                                                    y_proba_nb, y_proba_cox)

    km_results = {}

    if HAS_LIFELINES:
        log("\nKaplan-Meier Analysis:")

        for model_name, risk_col in [('NB', 'NB_Risk_Group'), ('Cox', 'Cox_Risk_Group')]:
            high_mask = df_risk[risk_col] == 'High'
            low_mask = df_risk[risk_col] == 'Low'

            if high_mask.sum() > 0 and low_mask.sum() > 0:
                kmf_high = KaplanMeierFitter()
                kmf_low = KaplanMeierFitter()

                kmf_high.fit(time[high_mask], event[high_mask], label='High Risk')
                kmf_low.fit(time[low_mask], event[low_mask], label='Low Risk')

                lr_result = logrank_test(time[high_mask], time[low_mask],
                                        event[high_mask], event[low_mask])

                median_high = kmf_high.median_survival_time_
                median_low = kmf_low.median_survival_time_

                # Calculate hazard ratio using Cox regression
                hr = np.nan
                hr_ci_low = np.nan
                hr_ci_high = np.nan
                try:
                    # Create dataframe for Cox regression
                    cox_df = pd.DataFrame({
                        'time': np.concatenate([time[high_mask], time[low_mask]]),
                        'event': np.concatenate([event[high_mask], event[low_mask]]),
                        'high_risk': np.concatenate([np.ones(high_mask.sum()), np.zeros(low_mask.sum())])
                    })
                    cph = CoxPHFitter()
                    cph.fit(cox_df, duration_col='time', event_col='event')
                    hr = np.exp(cph.params_['high_risk'])
                    hr_ci_low = np.exp(cph.confidence_intervals_.loc['high_risk', '95% lower-bound'])
                    hr_ci_high = np.exp(cph.confidence_intervals_.loc['high_risk', '95% upper-bound'])
                except Exception as e:
                    log(f"    Warning: Could not compute HR for {model_name}: {e}")

                km_results[model_name] = {
                    'n_high': high_mask.sum(),
                    'n_low': low_mask.sum(),
                    'median_pfs_high': median_high,
                    'median_pfs_low': median_low,
                    'logrank_stat': lr_result.test_statistic,
                    'logrank_p': lr_result.p_value,
                    'kmf_high': kmf_high,
                    'kmf_low': kmf_low,
                    'hazard_ratio': hr,
                    'hr_ci_low': hr_ci_low,
                    'hr_ci_high': hr_ci_high
                }

                log(f"\n  {model_name} Model:")
                log(f"    High Risk: n={high_mask.sum()}, median PFS={median_high:.1f} days")
                log(f"    Low Risk: n={low_mask.sum()}, median PFS={median_low:.1f} days")
                log(f"    Log-rank p-value: {lr_result.p_value:.4f}")

    return {
        'risk_df': df_risk,
        'agreement': agreement,
        'kappa': kappa,
        'agreement_analysis': agreement_analysis,
        'km_results': km_results,
        'nb_classifier': nb_clf,
        'cox_classifier': cox_clf
    }


def characterize_disagreements(y_true: np.ndarray, y_pred_nb: np.ndarray,
                                y_pred_cox: np.ndarray, y_proba_nb: np.ndarray,
                                y_proba_cox: np.ndarray) -> Dict:
    """Analyze disagreements between NB and Cox models."""

    agree = y_pred_nb == y_pred_cox
    disagree = ~agree

    # Types of disagreements
    nb_pos_cox_neg = (y_pred_nb == 1) & (y_pred_cox == 0)
    nb_neg_cox_pos = (y_pred_nb == 0) & (y_pred_cox == 1)

    # Who is correct when they disagree?
    nb_correct_when_disagree = disagree & (y_pred_nb == y_true)
    cox_correct_when_disagree = disagree & (y_pred_cox == y_true)

    # Probability differences in disagreements
    prob_diff_disagree = np.abs(y_proba_nb[disagree] - y_proba_cox[disagree]) if disagree.sum() > 0 else []

    return {
        'n_agree': agree.sum(),
        'n_disagree': disagree.sum(),
        'pct_agree': agree.mean() * 100,
        'n_nb_pos_cox_neg': nb_pos_cox_neg.sum(),
        'n_nb_neg_cox_pos': nb_neg_cox_pos.sum(),
        'nb_correct_when_disagree': nb_correct_when_disagree.sum(),
        'cox_correct_when_disagree': cox_correct_when_disagree.sum(),
        'mean_prob_diff_disagree': np.mean(prob_diff_disagree) if len(prob_diff_disagree) > 0 else 0
    }

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def plot_comparison_bars(summary_nb: Dict, summary_cox: Dict, summary_zeror: Dict,
                         paired_tests_nb_vs_zeror: Dict, paired_tests_cox_vs_zeror: Dict,
                         cohens_d_nb_vs_zeror: Dict, cohens_d_cox_vs_zeror: Dict,
                         title: str, filename: str,
                         paired_tests_nb_vs_cox: Dict = None, cohens_d_nb_vs_cox: Dict = None):
    """Create publication-quality bar chart comparing NB, Cox, and ZeroR.

    Shows all three pairwise significance comparisons with Bonferroni correction.
    """

    metrics = ['Accuracy', 'Balanced_Accuracy', 'F1', 'Precision', 'Recall',
               'Specificity', 'NPV', 'ROC_AUC', 'PRC_AUC']

    fig, axes = plt.subplots(3, 3, figsize=(16, 16))  # Increased height for more y-space
    axes = axes.flatten()

    def get_stars_bonferroni(p_val, n_comparisons=3):
        """Get significance stars with Bonferroni correction."""
        # Bonferroni-adjusted thresholds
        adj_001 = 0.001 / n_comparisons
        adj_01 = 0.01 / n_comparisons
        adj_05 = 0.05 / n_comparisons

        if p_val < adj_001:
            return '***'
        elif p_val < adj_01:
            return '**'
        elif p_val < adj_05:
            return '*'
        else:
            return 'ns'

    for idx, metric in enumerate(metrics):
        ax = axes[idx]

        nb_val, nb_lo, nb_hi = summary_nb.get(metric, (0.5, 0.5, 0.5))
        cox_val, cox_lo, cox_hi = summary_cox.get(metric, (0.5, 0.5, 0.5))
        zeror_val, zeror_lo, zeror_hi = summary_zeror.get(metric, (0.5, 0.5, 0.5))

        x = [0, 1, 2]
        values = [nb_val, cox_val, zeror_val]
        errors_low = [nb_val - nb_lo, cox_val - cox_lo, zeror_val - zeror_lo]
        errors_high = [nb_hi - nb_val, cox_hi - cox_val, zeror_hi - zeror_val]
        colors = [COLORMAP['NB'], COLORMAP['Cox'], COLORMAP['ZeroR']]

        bars = ax.bar(x, values, color=colors, edgecolor='black', linewidth=0.5)

        ax.errorbar(x, values, yerr=[errors_low, errors_high],
                   fmt='none', color='black', capsize=5, capthick=1.5)

        # Calculate vertical positions with proper spacing
        max_hi = max(nb_hi, cox_hi, zeror_hi)

        # Position mean value text above the error bar
        for i, (val, lo, hi) in enumerate([(nb_val, nb_lo, nb_hi),
                                           (cox_val, cox_lo, cox_hi),
                                           (zeror_val, zeror_lo, zeror_hi)]):
            # Mean value - positioned just above the bar
            mean_y = val + 0.02
            ax.text(i, mean_y, f'{val:.3f}',
                   ha='center', va='bottom', fontsize=9, fontweight='bold',
                   color='darkblue' if i < 2 else 'darkred')

            # CI text - positioned above the error bar whisker with more spacing
            ci_y = hi + 0.06
            ax.text(i, ci_y, f'[{lo:.2f}-{hi:.2f}]',
                   ha='center', va='bottom', fontsize=7, color='gray')

        # Get p-values for all three comparisons
        p_nb_zeror = paired_tests_nb_vs_zeror.get(metric, {}).get('p_value', 1)
        p_cox_zeror = paired_tests_cox_vs_zeror.get(metric, {}).get('p_value', 1)
        p_nb_cox = paired_tests_nb_vs_cox.get(metric, {}).get('p_value', 1) if paired_tests_nb_vs_cox else 1

        # Get significance stars with Bonferroni correction (3 comparisons)
        stars_nb_zeror = get_stars_bonferroni(p_nb_zeror, 3)
        stars_cox_zeror = get_stars_bonferroni(p_cox_zeror, 3)
        stars_nb_cox = get_stars_bonferroni(p_nb_cox, 3)

        # Calculate bracket positions with proper vertical stacking
        bracket_base = max_hi + 0.18  # Start brackets well above CI text
        bracket_spacing = 0.09  # Vertical space between brackets

        # Draw significance brackets for all three comparisons
        # Bracket 1: NB vs Cox (positions 0 and 1) - lowest bracket
        if stars_nb_cox != 'ns':
            y1 = bracket_base
            ax.plot([0, 0, 1, 1], [y1-0.015, y1, y1, y1-0.015], 'k-', linewidth=0.8)
            ax.text(0.5, y1 + 0.005, stars_nb_cox, ha='center', va='bottom', fontsize=10, color='purple')

        # Bracket 2: Cox vs ZeroR (positions 1 and 2) - middle bracket
        if stars_cox_zeror != 'ns':
            y2 = bracket_base + bracket_spacing
            ax.plot([1, 1, 2, 2], [y2-0.015, y2, y2, y2-0.015], 'k-', linewidth=0.8)
            ax.text(1.5, y2 + 0.005, stars_cox_zeror, ha='center', va='bottom', fontsize=10, color='darkgreen')

        # Bracket 3: NB vs ZeroR (positions 0 and 2) - top bracket (spans widest)
        if stars_nb_zeror != 'ns':
            y3 = bracket_base + 2 * bracket_spacing
            ax.plot([0, 0, 2, 2], [y3-0.015, y3, y3, y3-0.015], 'k-', linewidth=0.8)
            ax.text(1, y3 + 0.005, stars_nb_zeror, ha='center', va='bottom', fontsize=10, color='darkblue')

        ax.set_ylabel('Score')
        ax.set_title(metric.replace('_', ' '), fontsize=11, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(['NB\n(ML)', 'Cox\n(Survival)', 'ZeroR\n(Baseline)'])
        ax.set_ylim(0, 1.55)  # Increased y-limit for bracket space
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    legend_elements = [
        mpatches.Patch(facecolor=COLORMAP['NB'], edgecolor='black', label='NB (ML Champion)'),
        mpatches.Patch(facecolor=COLORMAP['Cox'], edgecolor='black', label='Cox (Survival)'),
        mpatches.Patch(facecolor=COLORMAP['ZeroR'], edgecolor='black', label='ZeroR (Baseline)')
    ]
    fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98))

    plt.suptitle(f'{title}\nBonferroni-corrected: *** p<0.001/3, ** p<0.01/3, * p<0.05/3, ns: not significant',
                fontsize=13, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    savefig(filename)

    # Store underlying data
    comparison_data = []
    for metric in metrics:
        nb_val, nb_lo, nb_hi = summary_nb.get(metric, (0.5, 0.5, 0.5))
        cox_val, cox_lo, cox_hi = summary_cox.get(metric, (0.5, 0.5, 0.5))
        zeror_val, zeror_lo, zeror_hi = summary_zeror.get(metric, (0.5, 0.5, 0.5))
        comparison_data.append({
            'Metric': metric, 'NB_Mean': nb_val, 'NB_CI_Low': nb_lo, 'NB_CI_High': nb_hi,
            'Cox_Mean': cox_val, 'Cox_CI_Low': cox_lo, 'Cox_CI_High': cox_hi,
            'ZeroR_Mean': zeror_val, 'ZeroR_CI_Low': zeror_lo, 'ZeroR_CI_High': zeror_hi
        })


def plot_head_to_head_comparison(summary_nb: Dict, summary_cox: Dict,
                                  paired_tests: Dict, cohens_d: Dict,
                                  title: str, filename: str):
    """Direct NB vs Cox comparison with difference visualization."""

    metrics = ['Accuracy', 'Balanced_Accuracy', 'F1', 'Precision', 'Recall',
               'Specificity', 'NPV', 'ROC_AUC', 'PRC_AUC', 'MCC']

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    ax1 = axes[0]
    x = np.arange(len(metrics))
    width = 0.35

    nb_vals = [summary_nb.get(m, (0.5, 0.5, 0.5))[0] for m in metrics]
    cox_vals = [summary_cox.get(m, (0.5, 0.5, 0.5))[0] for m in metrics]

    nb_errs = [[summary_nb.get(m, (0.5, 0.5, 0.5))[0] - summary_nb.get(m, (0.5, 0.5, 0.5))[1] for m in metrics],
               [summary_nb.get(m, (0.5, 0.5, 0.5))[2] - summary_nb.get(m, (0.5, 0.5, 0.5))[0] for m in metrics]]
    cox_errs = [[summary_cox.get(m, (0.5, 0.5, 0.5))[0] - summary_cox.get(m, (0.5, 0.5, 0.5))[1] for m in metrics],
                [summary_cox.get(m, (0.5, 0.5, 0.5))[2] - summary_cox.get(m, (0.5, 0.5, 0.5))[0] for m in metrics]]

    bars1 = ax1.bar(x - width/2, nb_vals, width, label='NB (ML)', color=COLORMAP['NB'],
                   yerr=nb_errs, capsize=3)
    bars2 = ax1.bar(x + width/2, cox_vals, width, label='Cox', color=COLORMAP['Cox'],
                   yerr=cox_errs, capsize=3)

    for i, metric in enumerate(metrics):
        p_val = paired_tests.get(metric, {}).get('p_value', 1)
        marker = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else ''

        if marker:
            # Position above error bars, not just bar height
            nb_top = nb_vals[i] + nb_errs[1][i]  # bar + upper error
            cox_top = cox_vals[i] + cox_errs[1][i]
            y_pos = max(nb_top, cox_top) + 0.03
            ax1.text(i, y_pos, marker, ha='center', fontsize=12, fontweight='bold')

    ax1.set_ylabel('Score')
    ax1.set_title('Performance Comparison: NB vs Cox')
    ax1.set_xticks(x)
    ax1.set_xticklabels([m.replace('_', '\n') for m in metrics], rotation=45, ha='right')
    ax1.legend()
    ax1.set_ylim(0, 1.2)
    ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

    ax2 = axes[1]
    differences = [nb_vals[i] - cox_vals[i] for i in range(len(metrics))]
    colors = [COLORMAP['NB'] if d > 0 else COLORMAP['Cox'] for d in differences]

    bars = ax2.barh(metrics, differences, color=colors, edgecolor='black')
    ax2.axvline(x=0, color='black', linewidth=1)

    for i, metric in enumerate(metrics):
        d = cohens_d.get(metric, 0)
        ax2.text(differences[i] + 0.01 if differences[i] >= 0 else differences[i] - 0.01,
                i, f'd={d:.2f}', va='center',
                ha='left' if differences[i] >= 0 else 'right',
                fontsize=9, style='italic')

    ax2.set_xlabel('Difference (NB - Cox)')
    ax2.set_title('Performance Difference\n(Positive = NB better)')
    ax2.set_xlim(-0.3, 0.3)

    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    savefig(filename)

    # Store underlying data
    h2h_data = []
    for i, metric in enumerate(metrics):
        h2h_data.append({
            'Metric': metric, 'NB_Value': nb_vals[i], 'Cox_Value': cox_vals[i],
            'Difference': differences[i], 'Cohen_d': cohens_d.get(metric, 0),
            'P_Value': paired_tests.get(metric, {}).get('p_value', np.nan)
        })


def plot_bootstrap_distributions(df_nb: pd.DataFrame, df_cox: pd.DataFrame,
                                  df_zeror: pd.DataFrame, filename: str):
    """Plot violin plots of bootstrap distributions."""

    metrics = ['F1', 'ROC_AUC', 'Balanced_Accuracy', 'MCC']

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for idx, metric in enumerate(metrics):
        ax = axes[idx]

        data = []
        for val in df_nb[metric]:
            data.append({'Model': 'NB (ML)', 'Value': val})
        for val in df_cox[metric]:
            data.append({'Model': 'Cox', 'Value': val})
        for val in df_zeror[metric]:
            data.append({'Model': 'ZeroR', 'Value': val})

        df_plot = pd.DataFrame(data)

        colors = [COLORMAP['NB'], COLORMAP['Cox'], COLORMAP['ZeroR']]

        parts = ax.violinplot([df_nb[metric].values, df_cox[metric].values, df_zeror[metric].values],
                              positions=[0, 1, 2], showmeans=True, showmedians=True)

        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(colors[i])
            pc.set_alpha(0.7)

        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['NB (ML)', 'Cox', 'ZeroR'])
        ax.set_ylabel(metric)
        ax.set_title(f'{metric} Bootstrap Distribution')
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

    plt.suptitle('Bootstrap Distributions Across Models', fontsize=14, fontweight='bold')
    plt.tight_layout()
    savefig(filename)

    # Store underlying data - full bootstrap distributions
    boot_data = []
    for metric in metrics:
        for i in range(len(df_nb)):
            boot_data.append({'Metric': metric, 'Model': 'NB', 'Bootstrap_Idx': i, 'Value': df_nb[metric].iloc[i]})
            boot_data.append({'Metric': metric, 'Model': 'Cox', 'Bootstrap_Idx': i, 'Value': df_cox[metric].iloc[i]})
            boot_data.append({'Metric': metric, 'Model': 'ZeroR', 'Bootstrap_Idx': i, 'Value': df_zeror[metric].iloc[i]})


def plot_confusion_matrices(y_true: np.ndarray, y_pred_nb: np.ndarray,
                            y_pred_cox: np.ndarray, filename: str):
    """Plot confusion matrices for NB and Cox with raw counts and row-normalized percentages."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for idx, (y_pred, title, color) in enumerate([
        (y_pred_nb, 'NB (ML) Champion', COLORMAP['NB']),
        (y_pred_cox, 'Cox Model', COLORMAP['Cox'])
    ]):
        ax = axes[idx]

        # Compute confusion matrix
        cm = confusion_matrix(y_true, y_pred, labels=[0, 1])

        # Row-normalize (each row sums to 100%)
        cm_row_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True) * 100

        # Create annotation labels with count and percentage
        annot_labels = np.empty_like(cm, dtype=object)
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                annot_labels[i, j] = f'{cm[i, j]}\n({cm_row_norm[i, j]:.1f}%)'

        # Plot heatmap with row-normalized values for color intensity
        sns.heatmap(cm_row_norm, annot=annot_labels, fmt='', cmap='Blues', ax=ax,
                   xticklabels=['Long PFS (L)', 'Short PFS (S)'],
                   yticklabels=['Long PFS (L)', 'Short PFS (S)'],
                   annot_kws={'size': 14, 'fontweight': 'bold'},
                   vmin=0, vmax=100, cbar_kws={'label': 'Row %'})

        ax.set_xlabel('Predicted', fontsize=12)
        ax.set_ylabel('Actual', fontsize=12)
        ax.set_title(f'{title}\n(count / row %)', fontsize=12)

        # Add row totals
        row_totals = cm.sum(axis=1)
        ax.text(2.3, 0.5, f'n={row_totals[0]}', va='center', fontsize=10, style='italic')
        ax.text(2.3, 1.5, f'n={row_totals[1]}', va='center', fontsize=10, style='italic')

    plt.suptitle('Confusion Matrix Comparison\n(Raw Count and Row-Normalized %)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    savefig(filename)

    # Store underlying data with percentages
    cm_nb = confusion_matrix(y_true, y_pred_nb, labels=[0, 1])
    cm_cox = confusion_matrix(y_true, y_pred_cox, labels=[0, 1])
    cm_nb_pct = cm_nb.astype(float) / cm_nb.sum(axis=1, keepdims=True) * 100
    cm_cox_pct = cm_cox.astype(float) / cm_cox.sum(axis=1, keepdims=True) * 100

    cm_data = {
        'NB_TN': cm_nb[0, 0], 'NB_FP': cm_nb[0, 1], 'NB_FN': cm_nb[1, 0], 'NB_TP': cm_nb[1, 1],
        'NB_TN_pct': cm_nb_pct[0, 0], 'NB_FP_pct': cm_nb_pct[0, 1],
        'NB_FN_pct': cm_nb_pct[1, 0], 'NB_TP_pct': cm_nb_pct[1, 1],
        'Cox_TN': cm_cox[0, 0], 'Cox_FP': cm_cox[0, 1], 'Cox_FN': cm_cox[1, 0], 'Cox_TP': cm_cox[1, 1],
        'Cox_TN_pct': cm_cox_pct[0, 0], 'Cox_FP_pct': cm_cox_pct[0, 1],
        'Cox_FN_pct': cm_cox_pct[1, 0], 'Cox_TP_pct': cm_cox_pct[1, 1]
    }


def plot_decision_curves(y_true: np.ndarray, y_proba_nb: np.ndarray,
                         y_proba_cox: np.ndarray, filename: str):
    """Plot decision curve analysis."""

    thresholds = np.linspace(0.01, 0.99, 100)

    dca_nb = decision_curve_analysis(y_true, y_proba_nb, thresholds)
    dca_cox = decision_curve_analysis(y_true, y_proba_cox, thresholds)

    fig, ax = plt.subplots(figsize=(10, 8))

    ax.plot(dca_nb['threshold'], dca_nb['net_benefit'],
           color=COLORMAP['NB'], linewidth=2, label='NB (ML)')
    ax.plot(dca_cox['threshold'], dca_cox['net_benefit'],
           color=COLORMAP['Cox'], linewidth=2, label='Cox')
    ax.plot(dca_nb['threshold'], dca_nb['treat_all'],
           color='gray', linewidth=1, linestyle='--', label='Treat All')
    ax.axhline(y=0, color='black', linewidth=1, linestyle='-', label='Treat None')

    ax.set_xlabel('Threshold Probability')
    ax.set_ylabel('Net Benefit')
    ax.set_title('Decision Curve Analysis')
    ax.legend(loc='upper right')
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.1, 0.6)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    savefig(filename)

    # Store underlying data
    dca_data = dca_nb[['threshold', 'net_benefit', 'treat_all']].copy()
    dca_data.columns = ['Threshold', 'NB_Net_Benefit', 'Treat_All']
    dca_data['Cox_Net_Benefit'] = dca_cox['net_benefit'].values


def plot_calibration_curves(y_true: np.ndarray, y_proba_nb: np.ndarray,
                            y_proba_cox: np.ndarray, filename: str):
    """Plot calibration curves for NB and Cox models."""

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot([0, 1], [0, 1], 'k--', label='Perfect calibration')

    prob_true_nb, prob_pred_nb = calibration_curve(y_true, y_proba_nb, n_bins=10)
    ax.plot(prob_pred_nb, prob_true_nb, 's-', color=COLORMAP['NB'],
           label='NB (ML)', linewidth=2, markersize=8)

    prob_true_cox, prob_pred_cox = calibration_curve(y_true, y_proba_cox, n_bins=10)
    ax.plot(prob_pred_cox, prob_true_cox, 'o-', color=COLORMAP['Cox'],
           label='Cox', linewidth=2, markersize=8)

    ax.set_xlabel('Mean Predicted Probability')
    ax.set_ylabel('Fraction of Positives')
    ax.set_title('Calibration Curves')
    ax.legend(loc='lower right')
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    savefig(filename)

    # Store underlying data (handle different array lengths)
    calib_rows = []
    for i, (pred, actual) in enumerate(zip(prob_pred_nb, prob_true_nb)):
        calib_rows.append({'Model': 'NB', 'Bin': i, 'Predicted': pred, 'Actual': actual})
    for i, (pred, actual) in enumerate(zip(prob_pred_cox, prob_true_cox)):
        calib_rows.append({'Model': 'Cox', 'Bin': i, 'Predicted': pred, 'Actual': actual})


def plot_roc_comparison(y_true: np.ndarray, y_proba_nb: np.ndarray,
                        y_proba_cox: np.ndarray, delong_result: Dict, filename: str):
    """Plot ROC curves with AUC comparison and DeLong test result."""

    fig, ax = plt.subplots(figsize=(8, 8))

    fpr_nb, tpr_nb, _ = roc_curve(y_true, y_proba_nb)
    auc_nb = roc_auc_score(y_true, y_proba_nb)
    ax.plot(fpr_nb, tpr_nb, color=COLORMAP['NB'], linewidth=2,
           label=f'NB (AUC={auc_nb:.3f})')

    fpr_cox, tpr_cox, _ = roc_curve(y_true, y_proba_cox)
    auc_cox = roc_auc_score(y_true, y_proba_cox)
    ax.plot(fpr_cox, tpr_cox, color=COLORMAP['Cox'], linewidth=2,
           label=f'Cox (AUC={auc_cox:.3f})')

    ax.plot([0, 1], [0, 1], 'k--', label='Random (AUC=0.5)')

    # Add DeLong test result
    delong_text = (f"DeLong Test:\n"
                   f"dAUC = {delong_result['difference']:.3f}\n"
                   f"z = {delong_result['z_score']:.2f}\n"
                   f"p = {delong_result['p_value']:.4f}")

    ax.text(0.55, 0.15, delong_text, transform=ax.transAxes,
           fontsize=10, verticalalignment='bottom',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlabel('False Positive Rate (1 - Specificity)')
    ax.set_ylabel('True Positive Rate (Sensitivity)')
    ax.set_title('ROC Curve Comparison with DeLong Test')
    ax.legend(loc='lower right')
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    savefig(filename)

    # Store underlying data
    roc_data = pd.DataFrame({
        'NB_FPR': fpr_nb, 'NB_TPR': tpr_nb,
        'Cox_FPR': np.interp(fpr_nb, fpr_cox, np.linspace(0, 1, len(fpr_cox))),
        'Cox_TPR': np.interp(fpr_nb, fpr_cox, tpr_cox)
    })
    roc_data['NB_AUC'] = auc_nb
    roc_data['Cox_AUC'] = auc_cox


def plot_nri_idi(nri_idi_results: Dict, filename: str):
    """Visualize NRI and IDI results with values displayed in bars."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: NRI components
    ax1 = axes[0]
    components = ['Events\n(up is good)', 'Non-Events\n(down is good)', 'Total NRI']
    values = [nri_idi_results['NRI_events'],
              nri_idi_results['NRI_non_events'],
              nri_idi_results['NRI']]
    colors = ['green' if v > 0 else 'red' for v in values]

    bars = ax1.bar(components, values, color=colors, edgecolor='black', alpha=0.8)
    ax1.axhline(y=0, color='black', linewidth=1)
    ax1.set_ylabel('NRI Component')
    ax1.set_title(f"Net Reclassification Improvement\nNRI = {nri_idi_results['NRI']:.3f} (p = {nri_idi_results['NRI_p']:.4f})")

    # Add values inside or above bars with better visibility
    for i, (bar, v) in enumerate(zip(bars, values)):
        # Position text inside bar if bar is tall enough, otherwise above
        bar_height = bar.get_height()
        if abs(bar_height) > 0.05:
            # Inside the bar
            y_pos = bar_height / 2
            text_color = 'white'
        else:
            # Above/below the bar
            y_pos = bar_height + 0.02 if bar_height >= 0 else bar_height - 0.02
            text_color = 'black'

        ax1.text(bar.get_x() + bar.get_width()/2, y_pos, f'{v:.3f}',
                ha='center', va='center', fontsize=12, fontweight='bold',
                color=text_color,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7) if abs(bar_height) <= 0.05 else None)

    # Set y-limits to ensure all values are visible
    y_min, y_max = ax1.get_ylim()
    ax1.set_ylim(min(y_min, -0.15), max(y_max, 0.15))

    # Right: IDI (Discrimination Slopes)
    ax2 = axes[1]
    categories = ['Old Model\n(Cox)', 'New Model\n(NB)']
    slopes = [nri_idi_results['slope_old'], nri_idi_results['slope_new']]

    bars = ax2.bar(categories, slopes, color=[COLORMAP['Cox'], COLORMAP['NB']], edgecolor='black', alpha=0.8)
    ax2.set_ylabel('Discrimination Slope')
    ax2.set_title(f"Integrated Discrimination Improvement\nIDI = {nri_idi_results['IDI']:.3f} (p = {nri_idi_results['IDI_p']:.4f})")

    # Add values inside bars
    for bar, v in zip(bars, slopes):
        bar_height = bar.get_height()
        # Position text inside the bar at 50% height
        y_pos = bar_height / 2
        ax2.text(bar.get_x() + bar.get_width()/2, y_pos, f'{v:.3f}',
                ha='center', va='center', fontsize=14, fontweight='bold',
                color='white')

    plt.suptitle('Reclassification Analysis: NB vs Cox', fontsize=14, fontweight='bold')
    plt.tight_layout()
    savefig(filename)

    # Store underlying data


def plot_km_comparison(km_results: Dict, filename: str):
    """Plot Kaplan-Meier curves for NB and Cox risk stratification with median droplines."""

    if not HAS_LIFELINES or not km_results:
        log("Skipping KM plot - lifelines not available or no results")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    for idx, (model_name, results) in enumerate(km_results.items()):
        ax = axes[idx]

        kmf_high = results['kmf_high']
        kmf_low = results['kmf_low']
        median_high = results['median_pfs_high']
        median_low = results['median_pfs_low']

        # Plot survival curves manually to control legend properly
        # Low Risk (green)
        sf_low = kmf_low.survival_function_
        ci_low = kmf_low.confidence_interval_survival_function_
        ax.step(sf_low.index, sf_low.iloc[:, 0], where='post', color='green', linewidth=2.5, label='Low Risk')
        ax.fill_between(ci_low.index, ci_low.iloc[:, 0], ci_low.iloc[:, 1],
                       step='post', alpha=0.25, color='green')

        # High Risk (red)
        sf_high = kmf_high.survival_function_
        ci_high = kmf_high.confidence_interval_survival_function_
        ax.step(sf_high.index, sf_high.iloc[:, 0], where='post', color='red', linewidth=2.5, label='High Risk')
        ax.fill_between(ci_high.index, ci_high.iloc[:, 0], ci_high.iloc[:, 1],
                       step='post', alpha=0.25, color='red')

        # Horizontal median line at y=0.5
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, linewidth=1)

        # Median droplines for Low Risk (green) - using KM estimate
        if not np.isnan(median_low) and not np.isinf(median_low):
            # Vertical line from median time down to x-axis
            ax.plot([median_low, median_low], [0, 0.5], color='green', linestyle=':', linewidth=1.5, alpha=0.8)
            # Horizontal line from y-axis to median time
            ax.plot([0, median_low], [0.5, 0.5], color='green', linestyle=':', linewidth=1.5, alpha=0.8)
            # Label the median
            ax.annotate(f'{median_low:.0f}d', xy=(median_low, 0.02), ha='center', va='bottom',
                       fontsize=9, color='darkgreen', fontweight='bold')

        # Median droplines for High Risk (red) - using KM estimate
        if not np.isnan(median_high) and not np.isinf(median_high):
            # Vertical line from median time down to x-axis
            ax.plot([median_high, median_high], [0, 0.5], color='red', linestyle=':', linewidth=1.5, alpha=0.8)
            # Horizontal line from y-axis to median time
            ax.plot([0, median_high], [0.5, 0.5], color='red', linestyle=':', linewidth=1.5, alpha=0.8)
            # Label the median
            ax.annotate(f'{median_high:.0f}d', xy=(median_high, 0.02), ha='center', va='bottom',
                       fontsize=9, color='darkred', fontweight='bold')

        # Stats text box (top right)
        stats_text = (f"High Risk (n={results['n_high']}): Median={median_high:.0f}d\n"
                     f"Low Risk (n={results['n_low']}): Median={median_low:.0f}d\n"
                     f"Log-rank p={results['logrank_p']:.4f}")

        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes,
               verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.9),
               fontsize=9, family='monospace')

        # Hazard ratio text box (bottom right, below legend area)
        hr = results.get('hazard_ratio', np.nan)
        hr_ci_low = results.get('hr_ci_low', np.nan)
        hr_ci_high = results.get('hr_ci_high', np.nan)

        if not np.isnan(hr):
            hr_text = f"HR = {hr:.2f} (95% CI: {hr_ci_low:.2f}-{hr_ci_high:.2f})"
            ax.text(0.5, -0.12, hr_text, transform=ax.transAxes,
                   ha='center', va='top', fontsize=11, fontweight='bold',
                   bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='orange', alpha=0.9))

        ax.set_title(f'{model_name} Model Risk Stratification', fontsize=12, fontweight='bold')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Survival Probability')

        # Legend with correct colors
        legend = ax.legend(loc='lower right', framealpha=0.9,
                          facecolor='white', edgecolor='gray')
        for text, color in zip(legend.get_texts(), ['green', 'red']):
            text.set_color(color)
            text.set_fontweight('bold')

        ax.set_ylim(0, 1.05)
        ax.set_xlim(0, None)

    plt.suptitle('Kaplan-Meier Survival by Risk Group', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])  # Leave room for HR text
    savefig(filename)

    # Store underlying data
    km_data = []
    for model_name, results in km_results.items():
        km_data.append({
            'Model': model_name, 'N_High': results['n_high'], 'N_Low': results['n_low'],
            'Median_PFS_High': results['median_pfs_high'], 'Median_PFS_Low': results['median_pfs_low'],
            'LogRank_Stat': results['logrank_stat'], 'LogRank_P': results['logrank_p'],
            'Hazard_Ratio': results.get('hazard_ratio', np.nan),
            'HR_CI_Low': results.get('hr_ci_low', np.nan),
            'HR_CI_High': results.get('hr_ci_high', np.nan)
        })


# =============================================================================
# NEW VISUALIZATIONS: ROC/PRC WITH CI, PATIENT STABILITY, FEATURE ANALYSIS
# =============================================================================

def plot_roc_with_confidence_bands(roc_curves: Dict, y_true: np.ndarray,
                                    y_proba_nb: np.ndarray, y_proba_cox: np.ndarray,
                                    delong_result: Dict, filename: str):
    """Plot ROC curves with bootstrap confidence bands (matching original style)."""

    fig, ax = plt.subplots(figsize=(10, 10))

    # Mean FPR axis for interpolation
    mean_fpr = np.linspace(0, 1, 100)

    # Process NB ROC curves
    tprs_nb = []
    for fpr, tpr in roc_curves.get('bootstrap_roc_nb', []):
        tpr_interp = np.interp(mean_fpr, fpr, tpr)
        tpr_interp[0] = 0.0
        tprs_nb.append(tpr_interp)

    if tprs_nb:
        tprs_nb = np.array(tprs_nb)
        mean_tpr_nb = np.mean(tprs_nb, axis=0)
        mean_tpr_nb[-1] = 1.0
        std_tpr_nb = np.std(tprs_nb, axis=0)
        tpr_upper_nb = np.clip(mean_tpr_nb + 1.96 * std_tpr_nb, 0, 1)
        tpr_lower_nb = np.clip(mean_tpr_nb - 1.96 * std_tpr_nb, 0, 1)

        auc_nb = roc_auc_score(y_true, y_proba_nb)
        ax.plot(mean_fpr, mean_tpr_nb, color=COLORMAP['NB'], linewidth=2.5,
               label=f'NB (ML) AUC={auc_nb:.3f}')
        ax.fill_between(mean_fpr, tpr_lower_nb, tpr_upper_nb,
                       color=COLORMAP['NB'], alpha=0.2, label='NB 95% CI')

    # Process Cox ROC curves
    tprs_cox = []
    for fpr, tpr in roc_curves.get('bootstrap_roc_cox', []):
        tpr_interp = np.interp(mean_fpr, fpr, tpr)
        tpr_interp[0] = 0.0
        tprs_cox.append(tpr_interp)

    if tprs_cox:
        tprs_cox = np.array(tprs_cox)
        mean_tpr_cox = np.mean(tprs_cox, axis=0)
        mean_tpr_cox[-1] = 1.0
        std_tpr_cox = np.std(tprs_cox, axis=0)
        tpr_upper_cox = np.clip(mean_tpr_cox + 1.96 * std_tpr_cox, 0, 1)
        tpr_lower_cox = np.clip(mean_tpr_cox - 1.96 * std_tpr_cox, 0, 1)

        auc_cox = roc_auc_score(y_true, y_proba_cox)
        ax.plot(mean_fpr, mean_tpr_cox, color=COLORMAP['Cox'], linewidth=2.5,
               label=f'Cox AUC={auc_cox:.3f}')
        ax.fill_between(mean_fpr, tpr_lower_cox, tpr_upper_cox,
                       color=COLORMAP['Cox'], alpha=0.2, label='Cox 95% CI')

    # Random line
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Random (AUC=0.5)')

    # DeLong test annotation
    delong_text = (f"DeLong Test\n"
                   f"dAUC={delong_result['difference']:.3f}\n"
                   f"z={delong_result['z_score']:.2f}\n"
                   f"p={delong_result['p_value']:.4f}")
    ax.text(0.55, 0.15, delong_text, transform=ax.transAxes, fontsize=11,
           verticalalignment='bottom',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))

    ax.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=12)
    ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=12)
    ax.set_title('ROC Curves with 95% Bootstrap Confidence Bands', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_aspect('equal')

    plt.tight_layout()
    savefig(filename)

    # Store underlying data for export


def plot_prc_with_confidence_bands(prc_curves: Dict, y_true: np.ndarray,
                                    y_proba_nb: np.ndarray, y_proba_cox: np.ndarray,
                                    filename: str):
    """Plot Precision-Recall curves with bootstrap confidence bands."""

    fig, ax = plt.subplots(figsize=(10, 10))

    # Mean recall axis for interpolation
    mean_recall = np.linspace(0, 1, 100)

    # Process NB PRC curves
    precs_nb = []
    for rec, prec in prc_curves.get('bootstrap_prc_nb', []):
        # Reverse for interpolation (recall is decreasing in original)
        prec_interp = np.interp(mean_recall, rec[::-1], prec[::-1])
        precs_nb.append(prec_interp)

    if precs_nb:
        precs_nb = np.array(precs_nb)
        mean_prec_nb = np.mean(precs_nb, axis=0)
        std_prec_nb = np.std(precs_nb, axis=0)
        prec_upper_nb = np.clip(mean_prec_nb + 1.96 * std_prec_nb, 0, 1)
        prec_lower_nb = np.clip(mean_prec_nb - 1.96 * std_prec_nb, 0, 1)

        ap_nb = average_precision_score(y_true, y_proba_nb)
        ax.plot(mean_recall, mean_prec_nb, color=COLORMAP['NB'], linewidth=2.5,
               label=f'NB (ML) AP={ap_nb:.3f}')
        ax.fill_between(mean_recall, prec_lower_nb, prec_upper_nb,
                       color=COLORMAP['NB'], alpha=0.2, label='NB 95% CI')

    # Process Cox PRC curves
    precs_cox = []
    for rec, prec in prc_curves.get('bootstrap_prc_cox', []):
        prec_interp = np.interp(mean_recall, rec[::-1], prec[::-1])
        precs_cox.append(prec_interp)

    if precs_cox:
        precs_cox = np.array(precs_cox)
        mean_prec_cox = np.mean(precs_cox, axis=0)
        std_prec_cox = np.std(precs_cox, axis=0)
        prec_upper_cox = np.clip(mean_prec_cox + 1.96 * std_prec_cox, 0, 1)
        prec_lower_cox = np.clip(mean_prec_cox - 1.96 * std_prec_cox, 0, 1)

        ap_cox = average_precision_score(y_true, y_proba_cox)
        ax.plot(mean_recall, mean_prec_cox, color=COLORMAP['Cox'], linewidth=2.5,
               label=f'Cox AP={ap_cox:.3f}')
        ax.fill_between(mean_recall, prec_lower_cox, prec_upper_cox,
                       color=COLORMAP['Cox'], alpha=0.2, label='Cox 95% CI')

    # Baseline (no-skill classifier)
    baseline = np.mean(y_true)
    ax.axhline(y=baseline, color='gray', linestyle='--', linewidth=1,
              label=f'No Skill (Prevalence={baseline:.2f})')

    ax.set_xlabel('Recall (Sensitivity)', fontsize=12)
    ax.set_ylabel('Precision (PPV)', fontsize=12)
    ax.set_title('Precision-Recall Curves with 95% Bootstrap Confidence Bands',
                fontsize=14, fontweight='bold')
    ax.legend(loc='lower left', fontsize=10)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

    plt.tight_layout()
    savefig(filename)

    # Store underlying data


def plot_patient_stability_analysis(patient_stability: Dict, filename: str):
    """Plot patient-level prediction stability analysis."""

    if not patient_stability:
        log("No patient stability data available")
        return

    df = pd.DataFrame.from_dict(patient_stability, orient='index')
    df['patient_id'] = df.index

    fig = plt.figure(figsize=(16, 12))

    # 1. Stability Distribution (NB vs Cox)
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.hist(df['stability_nb'], bins=20, alpha=0.7, color=COLORMAP['NB'], label='NB', edgecolor='black')
    ax1.hist(df['stability_cox'], bins=20, alpha=0.7, color=COLORMAP['Cox'], label='Cox', edgecolor='black')
    ax1.axvline(x=STABILITY_THRESHOLD, color='red', linestyle='--', linewidth=2, label=f'Threshold ({STABILITY_THRESHOLD})')
    ax1.set_xlabel('Prediction Stability')
    ax1.set_ylabel('Number of Patients')
    ax1.set_title('Prediction Stability Distribution')
    ax1.legend()

    # 2. Probability Standard Deviation
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.scatter(df['std_proba_nb'], df['std_proba_cox'], c=df['true_label'],
               cmap='RdYlBu', alpha=0.7, edgecolors='black', s=80)
    ax2.plot([0, 0.5], [0, 0.5], 'k--', alpha=0.5)
    ax2.set_xlabel('NB Probability Std Dev')
    ax2.set_ylabel('Cox Probability Std Dev')
    ax2.set_title('Prediction Uncertainty Comparison')
    ax2.set_xlim(0, 0.5)
    ax2.set_ylim(0, 0.5)

    # 3. Mean Probability vs True Label
    ax3 = fig.add_subplot(2, 3, 3)
    colors = ['green' if c else 'red' for c in df['correct_nb']]
    ax3.scatter(df['mean_proba_nb'], df['mean_proba_cox'], c=colors,
               alpha=0.7, edgecolors='black', s=80)
    ax3.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax3.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)
    ax3.set_xlabel('NB Mean Probability')
    ax3.set_ylabel('Cox Mean Probability')
    ax3.set_title('Mean Predictions (Green=Correct NB)')
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)

    # 4. Stability by True Label
    ax4 = fig.add_subplot(2, 3, 4)
    df_long = df[df['true_label'] == 0]
    df_short = df[df['true_label'] == 1]

    positions = [1, 2, 4, 5]
    data = [df_long['stability_nb'], df_short['stability_nb'],
            df_long['stability_cox'], df_short['stability_cox']]
    bp = ax4.boxplot(data, positions=positions, patch_artist=True, widths=0.6)

    colors = [COLORMAP['NB'], COLORMAP['NB'], COLORMAP['Cox'], COLORMAP['Cox']]
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax4.set_xticks([1.5, 4.5])
    ax4.set_xticklabels(['NB (ML)', 'Cox'])
    ax4.set_ylabel('Stability')
    ax4.set_title('Stability by Model and True Label\n(Left=Long PFS, Right=Short PFS)')
    ax4.axhline(y=STABILITY_THRESHOLD, color='red', linestyle='--', alpha=0.7)

    # 5. Agreement Heatmap (Square confusion matrix)
    ax5 = fig.add_subplot(2, 3, 5)
    agreement_matrix = np.zeros((2, 2))
    for _, row in df.iterrows():
        nb_pred = row['majority_vote_nb']
        cox_pred = row['majority_vote_cox']
        agreement_matrix[int(nb_pred), int(cox_pred)] += 1

    sns.heatmap(agreement_matrix.astype(int), annot=True, fmt='d', cmap='Blues', ax=ax5,
               xticklabels=['Long PFS', 'Short PFS'],
               yticklabels=['Long PFS', 'Short PFS'],
               annot_kws={'size': 14}, square=True, cbar_kws={'shrink': 0.8})
    ax5.set_xlabel('Cox Prediction')
    ax5.set_ylabel('NB Prediction')
    ax5.set_title('Model Agreement Matrix')
    ax5.set_aspect('equal')

    # 6. Summary Statistics
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.axis('off')

    n_stable_nb = df['is_stable_nb'].sum()
    n_stable_cox = df['is_stable_cox'].sum()
    n_correct_nb = df['correct_nb'].sum()
    n_correct_cox = df['correct_cox'].sum()
    n_total = len(df)

    summary_text = (
        f"PATIENT STABILITY SUMMARY\n"
        f"{'='*40}\n\n"
        f"Total Patients: {n_total}\n\n"
        f"NB (ML) Model:\n"
        f"  Stable predictions: {n_stable_nb}/{n_total} ({100*n_stable_nb/n_total:.1f}%)\n"
        f"  Correct (majority vote): {n_correct_nb}/{n_total} ({100*n_correct_nb/n_total:.1f}%)\n"
        f"  Mean stability: {df['stability_nb'].mean():.3f}\n\n"
        f"Cox Model:\n"
        f"  Stable predictions: {n_stable_cox}/{n_total} ({100*n_stable_cox/n_total:.1f}%)\n"
        f"  Correct (majority vote): {n_correct_cox}/{n_total} ({100*n_correct_cox/n_total:.1f}%)\n"
        f"  Mean stability: {df['stability_cox'].mean():.3f}\n\n"
        f"Model Agreement: {np.mean(df['majority_vote_nb'] == df['majority_vote_cox'])*100:.1f}%"
    )

    ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes,
            fontsize=11, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))

    plt.suptitle('Patient-Level Prediction Stability Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    savefig(filename)

    # Store underlying data

    return df


def plot_feature_correlation_analysis(X_nb: np.ndarray, X_cox: np.ndarray,
                                       y: np.ndarray, filename: str):
    """Plot feature correlation and importance analysis."""

    fig = plt.figure(figsize=(16, 12))

    # Combine all features
    all_features = NB_FEATURES + COX_FEATURES
    X_all = np.column_stack([X_nb, X_cox])
    df_features = pd.DataFrame(X_all, columns=all_features)

    # 1. Feature Correlation Heatmap
    ax1 = fig.add_subplot(2, 2, 1)
    if FEATURE_CORRELATION_METHOD == 'spearman':
        corr_matrix = df_features.corr(method='spearman')
    else:
        corr_matrix = df_features.corr(method='pearson')

    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
    sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.2f', cmap='RdBu_r',
               ax=ax1, vmin=-1, vmax=1, center=0,
               annot_kws={'size': 9})
    ax1.set_title(f'Feature Correlation ({FEATURE_CORRELATION_METHOD.capitalize()})')

    # 2. Feature Distributions by Class
    ax2 = fig.add_subplot(2, 2, 2)
    feature_means = []
    for i, feat in enumerate(all_features):
        class_0_mean = np.mean(X_all[y == 0, i])
        class_1_mean = np.mean(X_all[y == 1, i])
        feature_means.append({
            'Feature': feat,
            'Long PFS': class_0_mean,
            'Short PFS': class_1_mean,
            'Model': 'NB' if feat in NB_FEATURES else 'Cox'
        })

    df_means = pd.DataFrame(feature_means)

    x = np.arange(len(all_features))
    width = 0.35
    bars1 = ax2.bar(x - width/2, df_means['Long PFS'], width, label='Long PFS', color='green', alpha=0.7)
    bars2 = ax2.bar(x + width/2, df_means['Short PFS'], width, label='Short PFS', color='red', alpha=0.7)

    ax2.set_xticks(x)
    ax2.set_xticklabels(all_features, rotation=45, ha='right')
    ax2.set_ylabel('Mean Value (Standardized)')
    ax2.set_title('Feature Means by PFS Group')
    ax2.legend(loc='upper right')

    # Add model separation line
    ax2.axvline(x=len(NB_FEATURES) - 0.5, color='gray', linestyle='--', alpha=0.5)
    # Position model labels inside the plot (at 95% of y-max to avoid title overlap)
    y_label_pos = ax2.get_ylim()[1] * 0.92
    ax2.text(len(NB_FEATURES)/2 - 0.5, y_label_pos, 'NB Features',
            ha='center', fontsize=10, fontweight='bold', color=COLORMAP['NB'],
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=COLORMAP['NB'], alpha=0.8))
    ax2.text(len(NB_FEATURES) + len(COX_FEATURES)/2 - 0.5, y_label_pos, 'Cox Features',
            ha='center', fontsize=10, fontweight='bold', color=COLORMAP['Cox'],
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=COLORMAP['Cox'], alpha=0.8))

    # 3. Feature Effect Size (Cohen's d)
    ax3 = fig.add_subplot(2, 2, 3)
    cohens_d_features = []
    for i, feat in enumerate(all_features):
        x0 = X_all[y == 0, i]
        x1 = X_all[y == 1, i]
        pooled_std = np.sqrt(((len(x0)-1)*np.std(x0)**2 + (len(x1)-1)*np.std(x1)**2) / (len(x0)+len(x1)-2))
        d = (np.mean(x1) - np.mean(x0)) / pooled_std if pooled_std > 0 else 0
        cohens_d_features.append(d)

    colors = [COLORMAP['NB'] if f in NB_FEATURES else COLORMAP['Cox'] for f in all_features]
    bars = ax3.barh(all_features, cohens_d_features, color=colors, edgecolor='black', alpha=0.8)
    ax3.axvline(x=0, color='black', linewidth=1)
    ax3.axvline(x=0.8, color='red', linestyle='--', alpha=0.5, label='Large effect')
    ax3.axvline(x=-0.8, color='red', linestyle='--', alpha=0.5)
    ax3.set_xlabel("Cohen's d (Short - Long PFS)")
    ax3.set_title("Feature Effect Sizes\n(Positive = Higher in Short PFS)")
    ax3.legend()

    # 4. Cross-Model Feature Correlations
    ax4 = fig.add_subplot(2, 2, 4)

    # Calculate correlations between NB and Cox features
    cross_corr = np.zeros((len(NB_FEATURES), len(COX_FEATURES)))
    for i, nb_feat in enumerate(NB_FEATURES):
        for j, cox_feat in enumerate(COX_FEATURES):
            if FEATURE_CORRELATION_METHOD == 'spearman':
                corr, _ = stats.spearmanr(df_features[nb_feat], df_features[cox_feat])
            else:
                corr, _ = stats.pearsonr(df_features[nb_feat], df_features[cox_feat])
            cross_corr[i, j] = corr

    sns.heatmap(cross_corr, annot=True, fmt='.2f', cmap='RdBu_r',
               xticklabels=COX_FEATURES, yticklabels=NB_FEATURES,
               ax=ax4, vmin=-1, vmax=1, center=0,
               annot_kws={'size': 11})
    ax4.set_xlabel('Cox Features')
    ax4.set_ylabel('NB Features')
    ax4.set_title('Cross-Model Feature Correlations')

    plt.suptitle('Feature Correlation and Importance Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    savefig(filename)

    # Store underlying data
    feature_data = pd.DataFrame({
        'Feature': all_features,
        'Long_PFS_Mean': df_means['Long PFS'],
        'Short_PFS_Mean': df_means['Short PFS'],
        'Cohen_d': cohens_d_features,
        'Model': ['NB' if f in NB_FEATURES else 'Cox' for f in all_features]
    })


def plot_learning_curves(X_nb: np.ndarray, X_cox: np.ndarray, y: np.ndarray,
                          filename: str):
    """Plot learning curves for both models."""

    from sklearn.model_selection import learning_curve

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    train_sizes = LEARNING_CURVE_TRAIN_SIZES

    # NB Learning Curve
    ax1 = axes[0]
    try:
        nb_model = NBChampionClassifier(features=NB_FEATURES, calibrate=False)  # Uncalibrated for speed
        train_sizes_nb, train_scores_nb, val_scores_nb = learning_curve(
            nb_model, X_nb, y,
            train_sizes=train_sizes,
            cv=LEARNING_CURVE_CV_FOLDS,
            scoring='f1',
            n_jobs=-1,
            random_state=RANDOM_STATE
        )

        train_mean_nb = np.mean(train_scores_nb, axis=1)
        train_std_nb = np.std(train_scores_nb, axis=1)
        val_mean_nb = np.mean(val_scores_nb, axis=1)
        val_std_nb = np.std(val_scores_nb, axis=1)

        ax1.plot(train_sizes_nb, train_mean_nb, 'o-', color=COLORMAP['NB'],
                label='Training F1', linewidth=2)
        ax1.fill_between(train_sizes_nb, train_mean_nb - train_std_nb,
                        train_mean_nb + train_std_nb, alpha=0.2, color=COLORMAP['NB'])
        ax1.plot(train_sizes_nb, val_mean_nb, 's--', color='darkblue',
                label='Validation F1', linewidth=2)
        ax1.fill_between(train_sizes_nb, val_mean_nb - val_std_nb,
                        val_mean_nb + val_std_nb, alpha=0.2, color='darkblue')

        ax1.set_xlabel('Training Set Size')
        ax1.set_ylabel('F1 Score')
        ax1.set_title('NB (ML) Learning Curve')
        ax1.legend(loc='lower right')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, 1.1)
    except Exception as e:
        ax1.text(0.5, 0.5, f'Learning curve failed:\n{str(e)[:50]}',
                ha='center', va='center', transform=ax1.transAxes)

    # Cox Learning Curve
    ax2 = axes[1]
    try:
        cox_model = CoxRiskClassifier(betas=COX_BETAS, features=COX_FEATURES,
                                      threshold_method=RISK_THRESHOLD_METHOD)
        train_sizes_cox, train_scores_cox, val_scores_cox = learning_curve(
            cox_model, X_cox, y,
            train_sizes=train_sizes,
            cv=LEARNING_CURVE_CV_FOLDS,
            scoring='f1',
            n_jobs=-1,
            random_state=RANDOM_STATE
        )

        train_mean_cox = np.mean(train_scores_cox, axis=1)
        train_std_cox = np.std(train_scores_cox, axis=1)
        val_mean_cox = np.mean(val_scores_cox, axis=1)
        val_std_cox = np.std(val_scores_cox, axis=1)

        ax2.plot(train_sizes_cox, train_mean_cox, 'o-', color=COLORMAP['Cox'],
                label='Training F1', linewidth=2)
        ax2.fill_between(train_sizes_cox, train_mean_cox - train_std_cox,
                        train_mean_cox + train_std_cox, alpha=0.2, color=COLORMAP['Cox'])
        ax2.plot(train_sizes_cox, val_mean_cox, 's--', color='darkmagenta',
                label='Validation F1', linewidth=2)
        ax2.fill_between(train_sizes_cox, val_mean_cox - val_std_cox,
                        val_mean_cox + val_std_cox, alpha=0.2, color='darkmagenta')

        ax2.set_xlabel('Training Set Size')
        ax2.set_ylabel('F1 Score')
        ax2.set_title('Cox Learning Curve')
        ax2.legend(loc='lower right')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(0, 1.1)
    except Exception as e:
        ax2.text(0.5, 0.5, f'Learning curve failed:\n{str(e)[:50]}',
                ha='center', va='center', transform=ax2.transAxes)

    plt.suptitle('Learning Curves: Training Size vs Performance', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    savefig(filename)

    # Store underlying data
    learning_data = []
    if 'train_sizes_nb' in dir() and train_sizes_nb is not None:
        for i, size in enumerate(train_sizes_nb):
            learning_data.append({
                'Model': 'NB', 'Train_Size': size,
                'Train_F1_Mean': train_mean_nb[i] if 'train_mean_nb' in dir() else np.nan,
                'Train_F1_Std': train_std_nb[i] if 'train_std_nb' in dir() else np.nan,
                'Val_F1_Mean': val_mean_nb[i] if 'val_mean_nb' in dir() else np.nan,
                'Val_F1_Std': val_std_nb[i] if 'val_std_nb' in dir() else np.nan
            })
    if 'train_sizes_cox' in dir() and train_sizes_cox is not None:
        for i, size in enumerate(train_sizes_cox):
            learning_data.append({
                'Model': 'Cox', 'Train_Size': size,
                'Train_F1_Mean': train_mean_cox[i] if 'train_mean_cox' in dir() else np.nan,
                'Train_F1_Std': train_std_cox[i] if 'train_std_cox' in dir() else np.nan,
                'Val_F1_Mean': val_mean_cox[i] if 'val_mean_cox' in dir() else np.nan,
                'Val_F1_Std': val_std_cox[i] if 'val_std_cox' in dir() else np.nan
            })
    if learning_data:
        pass


def plot_comprehensive_metrics_comparison(summary_nb: Dict, summary_cox: Dict,
                                           summary_zeror: Dict, filename: str):
    """Plot comprehensive metrics comparison with all requested metrics."""

    # All metrics including extended ones
    metrics = ['Accuracy', 'Balanced_Accuracy', 'F1', 'Precision', 'Recall',
               'Specificity', 'NPV', 'ROC_AUC', 'PRC_AUC', 'MCC',
               'Kappa', 'Brier', 'ECE']

    available_metrics = [m for m in metrics if m in summary_nb]

    fig, ax = plt.subplots(figsize=(16, 8))

    x = np.arange(len(available_metrics))
    width = 0.25

    nb_vals = [summary_nb.get(m, (0.5, 0.5, 0.5))[0] for m in available_metrics]
    cox_vals = [summary_cox.get(m, (0.5, 0.5, 0.5))[0] for m in available_metrics]
    zeror_vals = [summary_zeror.get(m, (0.5, 0.5, 0.5))[0] for m in available_metrics]

    # Error bars
    nb_err_low = [summary_nb.get(m, (0.5, 0.5, 0.5))[0] - summary_nb.get(m, (0.5, 0.5, 0.5))[1]
                  for m in available_metrics]
    nb_err_high = [summary_nb.get(m, (0.5, 0.5, 0.5))[2] - summary_nb.get(m, (0.5, 0.5, 0.5))[0]
                   for m in available_metrics]
    cox_err_low = [summary_cox.get(m, (0.5, 0.5, 0.5))[0] - summary_cox.get(m, (0.5, 0.5, 0.5))[1]
                   for m in available_metrics]
    cox_err_high = [summary_cox.get(m, (0.5, 0.5, 0.5))[2] - summary_cox.get(m, (0.5, 0.5, 0.5))[0]
                    for m in available_metrics]

    bars1 = ax.bar(x - width, nb_vals, width, label='NB (ML)', color=COLORMAP['NB'],
                  yerr=[nb_err_low, nb_err_high], capsize=3, edgecolor='black')
    bars2 = ax.bar(x, cox_vals, width, label='Cox', color=COLORMAP['Cox'],
                  yerr=[cox_err_low, cox_err_high], capsize=3, edgecolor='black')
    bars3 = ax.bar(x + width, zeror_vals, width, label='ZeroR', color=COLORMAP['ZeroR'],
                  capsize=3, edgecolor='black', alpha=0.7)

    ax.set_ylabel('Score')
    ax.set_title('Comprehensive Metrics Comparison (Mean with 95% Percentile CI)')
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace('_', '\n') for m in available_metrics], rotation=45, ha='right')
    ax.legend(loc='upper right')
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylim(0, 1.1)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    savefig(filename)

    # Store underlying data
    comp_data = []
    for i, metric in enumerate(available_metrics):
        comp_data.append({
            'Metric': metric,
            'NB_Mean': nb_vals[i], 'NB_CI_Low': summary_nb.get(metric, (0.5, 0.5, 0.5))[1],
            'NB_CI_High': summary_nb.get(metric, (0.5, 0.5, 0.5))[2],
            'Cox_Mean': cox_vals[i], 'Cox_CI_Low': summary_cox.get(metric, (0.5, 0.5, 0.5))[1],
            'Cox_CI_High': summary_cox.get(metric, (0.5, 0.5, 0.5))[2],
            'ZeroR_Mean': zeror_vals[i]
        })

# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Main analysis pipeline."""

    log("="*70)
    log("COX vs ML CHAMPION COMPARISON ANALYSIS ")
    log("="*70)
    log(f"Output directory: {OUTPUT_DIR}")
    log(f"Bootstrap iterations: {N_BOOTSTRAP}")
    log(f"CV configuration: {CV_N_REPEATS}x{CV_N_FOLDS}-fold")
    log(f"Total evaluations: {TOTAL_EVALUATIONS}")

    # -------------------------------------------------------------------------
    # 1. LOAD DATA
    # -------------------------------------------------------------------------
    log("\n[1/8] Loading data...")

    if not os.path.exists(FULL_DATASET):
        log(f"[ERROR] Dataset not found: {FULL_DATASET}", "ERROR")
        log("Please update FULL_DATASET path in configuration.", "ERROR")
        return

    df = pd.read_excel(FULL_DATASET)
    log(f"  Loaded dataset: {df.shape}")

    X_nb = df[NB_FEATURES].values
    X_cox = df[COX_FEATURES].values
    y = df[LABEL_COLUMN].map(CLASS_MAP).values
    time = df[TIME_COLUMN].values
    event = df[EVENT_COLUMN].values
    patient_ids = df['patient_ID'].values if 'patient_ID' in df.columns else np.arange(len(y))

    log(f"  NB features: {NB_FEATURES}")
    log(f"  Cox features: {COX_FEATURES}")
    log(f"  Class distribution: {np.sum(y==0)} Long, {np.sum(y==1)} Short")

    # -------------------------------------------------------------------------
    # 2. BOOTSTRAP CV EVALUATION
    # -------------------------------------------------------------------------
    log("\n[2/9] Running Bootstrap CV evaluation...")

    cv_results = run_bootstrap_cv_evaluation(
        X_nb, X_cox, y,
        patient_ids=patient_ids,
        n_bootstrap=N_BOOTSTRAP,
        n_repeats=CV_N_REPEATS,
        n_folds=CV_N_FOLDS,
        seed=RANDOM_STATE
    )

    # -------------------------------------------------------------------------
    # 3. ADVANCED STATISTICAL TESTS
    # -------------------------------------------------------------------------
    log("\n[3/9] Running advanced statistical analyses...")

    preds = cv_results['all_predictions']

    # DeLong test
    log("  Running DeLong test for AUC comparison...")
    delong_result = delong_test(preds['y_true'], preds['y_proba_nb'], preds['y_proba_cox'])
    log(f"    DeLong: dAUC={delong_result['difference']:.3f}, p={delong_result['p_value']:.4f}")

    # NRI and IDI (NB vs Cox, treating Cox as reference)
    log("  Computing NRI and IDI...")
    nri_idi_result = compute_nri_idi(preds['y_true'], preds['y_proba_cox'],
                                     preds['y_proba_nb'], NRI_THRESHOLD)
    log(f"    NRI={nri_idi_result['NRI']:.3f} (p={nri_idi_result['NRI_p']:.4f})")
    log(f"    IDI={nri_idi_result['IDI']:.3f} (p={nri_idi_result['IDI_p']:.4f})")

    # Threshold optimization
    log("  Optimizing classification thresholds...")
    threshold_nb = optimize_threshold(preds['y_true'], preds['y_proba_nb'])
    threshold_cox = optimize_threshold(preds['y_true'], preds['y_proba_cox'])

    # Brier score decomposition
    log("  Computing Brier score decomposition...")
    brier_nb = brier_score_decomposition(preds['y_true'], preds['y_proba_nb'])
    brier_cox = brier_score_decomposition(preds['y_true'], preds['y_proba_cox'])

    # Debug logging for all advanced statistics

    # -------------------------------------------------------------------------
    # 4. VISUALIZATIONS
    # -------------------------------------------------------------------------
    log("\n[4/9] Generating visualizations...")

    # Main comparison bar chart (with all pairwise comparisons)
    plot_comparison_bars(
        cv_results['summary_nb'],
        cv_results['summary_cox'],
        cv_results['summary_zeror'],
        cv_results['paired_tests_nb_vs_zeror'],
        cv_results['paired_tests_cox_vs_zeror'],
        cv_results['cohens_d_nb_vs_zeror'],
        cv_results['cohens_d_cox_vs_zeror'],
        'Bootstrap CV: NB vs Cox vs ZeroR',
        'comparison_bars_all_models',
        paired_tests_nb_vs_cox=cv_results['paired_tests_nb_vs_cox'],
        cohens_d_nb_vs_cox=cv_results['cohens_d_nb_vs_cox']
    )

    # Head-to-head NB vs Cox
    plot_head_to_head_comparison(
        cv_results['summary_nb'],
        cv_results['summary_cox'],
        cv_results['paired_tests_nb_vs_cox'],
        cv_results['cohens_d_nb_vs_cox'],
        'Head-to-Head: NB (ML) vs Cox (Survival)',
        'head_to_head_nb_vs_cox'
    )

    # Bootstrap distributions
    plot_bootstrap_distributions(
        cv_results['bootstrap_metrics_nb'],
        cv_results['bootstrap_metrics_cox'],
        cv_results['bootstrap_metrics_zeror'],
        'bootstrap_distributions'
    )

    # Confusion matrices
    plot_confusion_matrices(
        preds['y_true'], preds['y_pred_nb'], preds['y_pred_cox'],
        'confusion_matrices'
    )

    # ROC comparison with DeLong
    plot_roc_comparison(
        preds['y_true'], preds['y_proba_nb'], preds['y_proba_cox'],
        delong_result, 'roc_comparison_delong'
    )

    # Calibration curves
    plot_calibration_curves(
        preds['y_true'], preds['y_proba_nb'], preds['y_proba_cox'],
        'calibration_curves'
    )

    # Decision curve analysis
    plot_decision_curves(
        preds['y_true'], preds['y_proba_nb'], preds['y_proba_cox'],
        'decision_curve_analysis'
    )

    # NRI/IDI visualization
    plot_nri_idi(nri_idi_result, 'nri_idi_analysis')

    # NEW: ROC with confidence bands
    log("  Generating ROC curve with confidence bands...")
    if cv_results.get('roc_curves'):
        plot_roc_with_confidence_bands(
            cv_results['roc_curves'],
            preds['y_true'], preds['y_proba_nb'], preds['y_proba_cox'],
            delong_result, 'roc_with_confidence_bands'
        )

    # NEW: PRC with confidence bands
    log("  Generating PRC curve with confidence bands...")
    if cv_results.get('prc_curves'):
        plot_prc_with_confidence_bands(
            cv_results['prc_curves'],
            preds['y_true'], preds['y_proba_nb'], preds['y_proba_cox'],
            'prc_with_confidence_bands'
        )

    # NEW: Patient stability analysis
    log("  Generating patient stability analysis...")
    patient_stability_df = None
    if cv_results.get('patient_tracking', {}).get('patient_stability'):
        patient_stability_df = plot_patient_stability_analysis(
            cv_results['patient_tracking']['patient_stability'],
            'patient_stability_analysis'
        )

    # NEW: Comprehensive metrics comparison
    log("  Generating comprehensive metrics comparison...")
    plot_comprehensive_metrics_comparison(
        cv_results['summary_nb'],
        cv_results['summary_cox'],
        cv_results['summary_zeror'],
        'comprehensive_metrics_comparison'
    )

    # NEW: Feature correlation analysis
    log("  Generating feature correlation analysis...")
    plot_feature_correlation_analysis(X_nb, X_cox, y, 'feature_correlation_analysis')

    # NEW: Learning curves
    log("  Generating learning curves...")
    plot_learning_curves(X_nb, X_cox, y, 'learning_curves')

    # -------------------------------------------------------------------------
    # 5. RISK STRATIFICATION
    # -------------------------------------------------------------------------
    log("\n[5/9] Running risk stratification analysis...")

    risk_results = run_risk_stratification(
        X_nb, X_cox, y, time, event, patient_ids
    )

    if risk_results['km_results']:
        plot_km_comparison(risk_results['km_results'], 'km_risk_stratification')

    # -------------------------------------------------------------------------
    # 6. PERMUTATION TESTS (Optional - computationally intensive)
    # -------------------------------------------------------------------------
    log("\n[6/9] Running permutation tests...")

    # NB permutation test
    perm_nb = run_permutation_test(
        X_nb, y,
        NBChampionClassifier,
        {'features': NB_FEATURES, 'calibrate': True},
        n_permutations=min(N_PERMUTATIONS, 500),  # Reduced for speed
        seed=RANDOM_STATE
    )

    # Cox permutation test
    perm_cox = run_permutation_test(
        X_cox, y,
        CoxRiskClassifier,
        {'betas': COX_BETAS, 'features': COX_FEATURES, 'threshold_method': RISK_THRESHOLD_METHOD},
        n_permutations=min(N_PERMUTATIONS, 500),
        seed=RANDOM_STATE
    )

    # -------------------------------------------------------------------------
    # 7. EXPORT RESULTS
    # -------------------------------------------------------------------------
    log("\n[7/9] Exporting results...")

    # Summary table
    summary_rows = []
    for metric in CORE_METRICS + ['MCC']:
        nb_val = cv_results['summary_nb'].get(metric, (np.nan, np.nan, np.nan))
        cox_val = cv_results['summary_cox'].get(metric, (np.nan, np.nan, np.nan))
        zeror_val = cv_results['summary_zeror'].get(metric, (np.nan, np.nan, np.nan))

        p_nb_zeror = cv_results['paired_tests_nb_vs_zeror'].get(metric, {}).get('p_value', np.nan)
        p_cox_zeror = cv_results['paired_tests_cox_vs_zeror'].get(metric, {}).get('p_value', np.nan)
        p_nb_cox = cv_results['paired_tests_nb_vs_cox'].get(metric, {}).get('p_value', np.nan)

        summary_rows.append({
            'Metric': metric,
            'NB_Mean': nb_val[0],
            'NB_CI_Lower': nb_val[1],
            'NB_CI_Upper': nb_val[2],
            'Cox_Mean': cox_val[0],
            'Cox_CI_Lower': cox_val[1],
            'Cox_CI_Upper': cox_val[2],
            'ZeroR_Mean': zeror_val[0],
            'ZeroR_CI_Lower': zeror_val[1],
            'ZeroR_CI_Upper': zeror_val[2],
            'P_NB_vs_ZeroR': p_nb_zeror,
            'P_Cox_vs_ZeroR': p_cox_zeror,
            'P_NB_vs_Cox': p_nb_cox,
            'Cohen_d_NB_vs_ZeroR': cv_results['cohens_d_nb_vs_zeror'].get(metric, np.nan),
            'Cohen_d_Cox_vs_ZeroR': cv_results['cohens_d_cox_vs_zeror'].get(metric, np.nan),
            'Cohen_d_NB_vs_Cox': cv_results['cohens_d_nb_vs_cox'].get(metric, np.nan)
        })

    df_summary = pd.DataFrame(summary_rows)

    # Advanced statistics summary
    df_advanced = pd.DataFrame([{
        'Analysis': 'DeLong Test',
        'Statistic': 'AUC Difference',
        'Value': delong_result['difference'],
        'Z_Score': delong_result['z_score'],
        'P_Value': delong_result['p_value'],
        'Significant': delong_result['significant']
    }, {
        'Analysis': 'NRI',
        'Statistic': 'Net Reclassification',
        'Value': nri_idi_result['NRI'],
        'Z_Score': nri_idi_result['NRI_z'],
        'P_Value': nri_idi_result['NRI_p'],
        'Significant': nri_idi_result['NRI_p'] < 0.05
    }, {
        'Analysis': 'IDI',
        'Statistic': 'Discrimination Improvement',
        'Value': nri_idi_result['IDI'],
        'Z_Score': nri_idi_result['IDI_z'],
        'P_Value': nri_idi_result['IDI_p'],
        'Significant': nri_idi_result['IDI_p'] < 0.05
    }])

    # Brier decomposition
    df_brier = pd.DataFrame([
        {'Model': 'NB', **brier_nb},
        {'Model': 'Cox', **brier_cox}
    ])

    # Threshold optimization
    df_threshold = pd.DataFrame([
        {'Model': 'NB', 'Method': 'Youden', **threshold_nb['youden']},
        {'Model': 'NB', 'Method': 'F1', **threshold_nb['f1']},
        {'Model': 'Cox', 'Method': 'Youden', **threshold_cox['youden']},
        {'Model': 'Cox', 'Method': 'F1', **threshold_cox['f1']}
    ])

    # Agreement analysis
    df_agreement = pd.DataFrame([risk_results['agreement_analysis']])

    # Patient stability data
    if patient_stability_df is not None:
        df_patient_stability = patient_stability_df
    elif cv_results.get('patient_tracking', {}).get('patient_stability'):
        df_patient_stability = pd.DataFrame.from_dict(
            cv_results['patient_tracking']['patient_stability'], orient='index'
        )
        df_patient_stability['patient_id'] = df_patient_stability.index
    else:
        df_patient_stability = pd.DataFrame()

    # DCA data export
    log("  Exporting Decision Curve Analysis data...")
    thresholds_dca = np.linspace(0.01, 0.99, 100)
    dca_nb_data = decision_curve_analysis(preds['y_true'], preds['y_proba_nb'], thresholds_dca)
    dca_cox_data = decision_curve_analysis(preds['y_true'], preds['y_proba_cox'], thresholds_dca)
    df_dca = pd.merge(
        dca_nb_data.rename(columns={'net_benefit': 'NB_net_benefit', 'treat_all': 'treat_all'}),
        dca_cox_data[['threshold', 'net_benefit']].rename(columns={'net_benefit': 'Cox_net_benefit'}),
        on='threshold'
    )

    # Calibration curve data
    log("  Exporting calibration curve data...")
    try:
        prob_true_nb, prob_pred_nb = calibration_curve(preds['y_true'], preds['y_proba_nb'], n_bins=10)
        prob_true_cox, prob_pred_cox = calibration_curve(preds['y_true'], preds['y_proba_cox'], n_bins=10)
        df_calibration = pd.DataFrame({
            'NB_predicted': prob_pred_nb,
            'NB_actual': prob_true_nb,
            'Cox_predicted': prob_pred_cox,
            'Cox_actual': prob_true_cox
        })
    except:
        df_calibration = pd.DataFrame()

    # Permutation test distributions
    log("  Exporting permutation test data...")
    perm_data_rows = []
    for metric in ['F1', 'ROC_AUC', 'Accuracy']:
        if metric in perm_nb.get('perm_distributions', {}):
            for i, val in enumerate(perm_nb['perm_distributions'][metric]):
                perm_data_rows.append({'model': 'NB', 'metric': metric, 'perm_idx': i, 'value': val})
        if metric in perm_cox.get('perm_distributions', {}):
            for i, val in enumerate(perm_cox['perm_distributions'][metric]):
                perm_data_rows.append({'model': 'Cox', 'metric': metric, 'perm_idx': i, 'value': val})
    df_perm_dist = pd.DataFrame(perm_data_rows) if perm_data_rows else pd.DataFrame()

    # All predictions raw data
    log("  Exporting raw predictions...")
    df_all_predictions = pd.DataFrame({
        'y_true': preds['y_true'],
        'y_pred_nb': preds['y_pred_nb'],
        'y_pred_cox': preds['y_pred_cox'],
        'y_proba_nb': preds['y_proba_nb'],
        'y_proba_cox': preds['y_proba_cox']
    })

    # Configuration summary
    log("  Creating configuration summary...")
    df_config = create_configuration_summary()

    # KM survival data (including hazard ratios)
    log("  Exporting KM curve data...")
    km_data_rows = []
    if risk_results['km_results']:
        for model, km in risk_results['km_results'].items():
            km_data_rows.append({
                'model': model,
                'group': 'High Risk',
                'n': km['n_high'],
                'median_pfs': km['median_pfs_high'],
                'logrank_stat': km['logrank_stat'],
                'logrank_p': km['logrank_p'],
                'hazard_ratio': km.get('hazard_ratio', np.nan),
                'hr_ci_low': km.get('hr_ci_low', np.nan),
                'hr_ci_high': km.get('hr_ci_high', np.nan)
            })
            km_data_rows.append({
                'model': model,
                'group': 'Low Risk',
                'n': km['n_low'],
                'median_pfs': km['median_pfs_low'],
                'logrank_stat': km['logrank_stat'],
                'logrank_p': km['logrank_p'],
                'hazard_ratio': km.get('hazard_ratio', np.nan),
                'hr_ci_low': km.get('hr_ci_low', np.nan),
                'hr_ci_high': km.get('hr_ci_high', np.nan)
            })
    df_km_summary = pd.DataFrame(km_data_rows) if km_data_rows else pd.DataFrame()

    # Export to comprehensive Excel
    log("  Creating comprehensive Excel file...")
    sheets = {
        '1_Configuration': df_config,
        '2_Summary_Metrics': df_summary,
        '3_Advanced_Stats': df_advanced,
        '4_Brier_Decomp': df_brier,
        '5_Threshold_Opt': df_threshold,
        '6_Model_Agreement': df_agreement,
        '7_NB_Bootstrap': cv_results['bootstrap_metrics_nb'],
        '8_Cox_Bootstrap': cv_results['bootstrap_metrics_cox'],
        '9_ZeroR_Bootstrap': cv_results['bootstrap_metrics_zeror'],
        '10_Risk_Strat': risk_results['risk_df'],
        '11_NRI_IDI': pd.DataFrame([nri_idi_result]),
        '12_Patient_Stability': df_patient_stability,
        '13_DCA_Data': df_dca,
        '14_Calibration_Data': df_calibration,
        '15_Raw_Predictions': df_all_predictions,
        '16_Perm_Distrib': df_perm_dist,
        '17_KM_Summary': df_km_summary,
    }

    save_excel(os.path.join(OUTPUT_DIR, 'COX_vs_ML_COMPREHENSIVE_RESULTS.xlsx'), sheets)

    # Also save separate files for large datasets
    log("  Saving separate data files...")

    # Raw predictions to CSV
    df_all_predictions.to_csv(os.path.join(OUTPUT_DIR, 'raw_predictions.csv'), index=False)

    # Bootstrap metrics to separate file
    bootstrap_sheets = {
        'NB_Bootstrap_Full': cv_results['bootstrap_metrics_nb'],
        'Cox_Bootstrap_Full': cv_results['bootstrap_metrics_cox'],
        'ZeroR_Bootstrap_Full': cv_results['bootstrap_metrics_zeror']
    }
    save_excel(os.path.join(OUTPUT_DIR, 'bootstrap_metrics_full.xlsx'), bootstrap_sheets)

    # -------------------------------------------------------------------------
    # 8. FINAL SUMMARY
    # -------------------------------------------------------------------------
    log("\n[8/9] Final Summary")
    log("="*70)

    print("\n" + "="*70)
    print("FINAL EVALUATION: COX vs ML Champion")
    print("="*70)

    print(f"\nModel Features:")
    print(f"  NB (ML):  {NB_FEATURES}")
    print(f"  Cox:      {COX_FEATURES}")
    print(f"  Overlap:  {set(NB_FEATURES) & set(COX_FEATURES) or 'None'}")

    print(f"\nKey Metrics (Mean [95% CI]):")
    for metric in ['F1', 'ROC_AUC', 'Balanced_Accuracy']:
        nb = cv_results['summary_nb'].get(metric, (np.nan, np.nan, np.nan))
        cox = cv_results['summary_cox'].get(metric, (np.nan, np.nan, np.nan))
        p = cv_results['paired_tests_nb_vs_cox'].get(metric, {}).get('p_value', 1)

        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        winner = 'NB' if nb[0] > cox[0] else 'Cox' if cox[0] > nb[0] else 'Tie'

        print(f"  {metric}:")
        print(f"    NB:  {nb[0]:.3f} [{nb[1]:.3f}-{nb[2]:.3f}]")
        print(f"    Cox: {cox[0]:.3f} [{cox[1]:.3f}-{cox[2]:.3f}]")
        print(f"    p={p:.4f} {sig} -> {winner}")

    print(f"\nAdvanced Comparisons:")
    print(f"  DeLong Test: dAUC={delong_result['difference']:.3f}, p={delong_result['p_value']:.4f}")
    print(f"  NRI: {nri_idi_result['NRI']:.3f}, p={nri_idi_result['NRI_p']:.4f}")
    print(f"  IDI: {nri_idi_result['IDI']:.3f}, p={nri_idi_result['IDI_p']:.4f}")

    print(f"\nRisk Stratification:")
    print(f"  Model agreement: {risk_results['agreement']:.1%}")
    print(f"  Cohen's Kappa:   {risk_results['kappa']:.3f}")

    if risk_results['km_results']:
        print(f"\nKaplan-Meier Validation:")
        for model, km in risk_results['km_results'].items():
            sig = '***' if km['logrank_p'] < 0.001 else '**' if km['logrank_p'] < 0.01 else '*' if km['logrank_p'] < 0.05 else 'ns'
            print(f"  {model}: Log-rank p={km['logrank_p']:.4f} {sig}")

    print(f"\nPermutation Test p-values:")
    print(f"  NB F1: p={perm_nb['p_values'].get('F1', np.nan):.4f}")
    print(f"  Cox F1: p={perm_cox['p_values'].get('F1', np.nan):.4f}")

    # Patient stability summary
    if not df_patient_stability.empty:
        print(f"\nPatient Prediction Stability:")
        n_total = len(df_patient_stability)
        n_stable_nb = df_patient_stability['is_stable_nb'].sum() if 'is_stable_nb' in df_patient_stability.columns else 0
        n_stable_cox = df_patient_stability['is_stable_cox'].sum() if 'is_stable_cox' in df_patient_stability.columns else 0
        print(f"  NB stable predictions: {n_stable_nb}/{n_total} ({100*n_stable_nb/n_total:.1f}%)")
        print(f"  Cox stable predictions: {n_stable_cox}/{n_total} ({100*n_stable_cox/n_total:.1f}%)")

    print("\n" + "="*70)
    print(f"Results saved to: {OUTPUT_DIR}")
    print("="*70)

    # -------------------------------------------------------------------------
    # 9. CLEANUP AND RETURN
    # -------------------------------------------------------------------------
    log("\n[9/9] Analysis complete!")
    log(f"Total visualizations generated: 15+")
    log(f"Excel file: COX_vs_ML_COMPREHENSIVE_RESULTS.xlsx")

    return {
        'cv_results': cv_results,
        'delong_result': delong_result,
        'nri_idi_result': nri_idi_result,
        'risk_results': risk_results,
        'perm_nb': perm_nb,
        'perm_cox': perm_cox,
        'patient_stability_df': df_patient_stability
    }


if __name__ == "__main__":
    main()
