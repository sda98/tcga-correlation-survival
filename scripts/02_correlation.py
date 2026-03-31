#!/usr/bin/env python3
"""
02_correlation.py

Computes Spearman correlation and orthogonal regression (TLS via SVD)
between two genes across TCGA pan-cancer and AML samples.
Produces scatter plots and saves high-expression samples.

"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib


# Use non-interactive backend 
matplotlib.use("Agg")


# ============================================================
# Configuration
# ============================================================

EXPRESSION_FILE = os.path.join("results", "expression_clean.tsv")
RESULTS_DIR = "results"

# Genes to analyze (Other genes can be used)
GENE1 = "ING5"
GENE2 = "FOXP1"

# AML samples start with this prefix
AML_PREFIX = "TCGA-AB"

# Plot axis limits
XLIM = (-0.1, 8.1)
YLIM = (-0.1, 8.1)


# ============================================================
# Helper functions
# ============================================================

def load_gene_pair(expression_file, gene1, gene2, sample_filter=None):
    """
    Extract two genes from the expression matrix and return a DataFrame
    with samples as rows and genes as columns.
    """
    # Load the full expression matrix
    tpm = pd.read_csv(expression_file, sep="\t", index_col=0)

    # Filter to the two genes
    expr = tpm.loc[[gene1, gene2]]

    # Filter to specific samples if requested
    if sample_filter:
        cols = [c for c in expr.columns if c.startswith(sample_filter)]
        expr = expr[cols]

    # Transpose so rows = samples, columns = genes
    expr_df = expr.T
    expr_df.columns = [gene1, gene2]

    # Drop any rows with NaN
    expr_df = expr_df.dropna()

    return expr_df


def spearman_correlation(x, y):
    """
    Compute Spearman rank correlation.
    Returns: (rho, p_value)
    """
    rho, pval = stats.spearmanr(x, y)
    return round(rho, 2), pval


def format_pvalue(pval):
    """
    Format p-value for display on plots.
    """
    if pval < 2.2e-16:
        return "< 2.2e-16"
    else:
        return f"= {pval:.2g}"


def tls_regression(x, y):
    """
    Total Least Squares (orthogonal) regression via SVD.
    Returns: (slope, intercept)
    """
    mx, my = np.mean(x), np.mean(y)
    M = np.column_stack([x - mx, y - my])
    _, _, Vt = np.linalg.svd(M, full_matrices=False)
    slope = Vt[0, 1] / Vt[0, 0]
    intercept = my - slope * mx
    return slope, intercept


def clip_line_to_box(slope, intercept, xlim, ylim):
    """
    Clip the regression line to fit within the plot box.
    Returns: (x_start, y_start, x_end, y_end)
    """
    if not np.isfinite(slope):
        # Vertical line case
        mx = -intercept / slope if np.isfinite(slope) else 0
        return (mx, ylim[0], mx, ylim[1])

    # Find where line intersects all four edges of the box
    candidates = []

    # Left edge: x = xlim[0]
    y_at_left = intercept + slope * xlim[0]
    candidates.append((xlim[0], y_at_left))

    # Right edge: x = xlim[1]
    y_at_right = intercept + slope * xlim[1]
    candidates.append((xlim[1], y_at_right))

    # Bottom edge: y = ylim[0]
    if slope != 0:
        x_at_bottom = (ylim[0] - intercept) / slope
        candidates.append((x_at_bottom, ylim[0]))

    # Top edge: y = ylim[1]
    if slope != 0:
        x_at_top = (ylim[1] - intercept) / slope
        candidates.append((x_at_top, ylim[1]))

    # Keep only points inside the box
    inside = [
        (cx, cy)
        for cx, cy in candidates
        if xlim[0] <= cx <= xlim[1] and ylim[0] <= cy <= ylim[1]
    ]

    if len(inside) < 2:
        return None

    # Inset slightly so line doesn't touch the frame
    p1 = np.array(inside[0])
    p2 = np.array(inside[1])
    direction = p2 - p1
    length = np.linalg.norm(direction)
    if length > 0:
        u = direction / length
        inset = 0.01 * length
        p1 = p1 + inset * u
        p2 = p2 - inset * u

    return (p1[0], p1[1], p2[0], p2[1])


def build_annotation(rho, pval_str):
    """
    Build the annotation string for the plot.
    """
    return f"ρ = {rho}\nP {pval_str}"


def make_scatter_plot(expr_df, gene1, gene2, slope, intercept, annotation,
                      xlim, ylim, point_color, title_suffix, output_path):
    """
    Create a scatter plot with TLS regression line and stats annotation.
    """
    fig, ax = plt.subplots(figsize=(8, 9.5))

    x = expr_df[gene1].values
    y = expr_df[gene2].values

    # Scatter points
    ax.scatter(x, y, s=50, facecolors=point_color, edgecolors="black",
               linewidths=0.6, zorder=2)

    # Regression line clipped to data range with padding
    x_min, x_max = float(x.min()), float(x.max())
    x_pad = 0.1 * (x_max - x_min)
    line_xlim = (x_min - x_pad, x_max + x_pad)
    line_coords = clip_line_to_box(slope, intercept, line_xlim, ylim)
    if line_coords:
        ax.plot(
            [line_coords[0], line_coords[2]],
            [line_coords[1], line_coords[3]],
            color="black", linewidth=3, solid_capstyle="round", zorder=3,
        )

    # Annotation text
    ax.text(xlim[0] + 0.2, ylim[1] - 0.5, annotation,
            fontsize=19, verticalalignment="top", horizontalalignment="left")
            
    # Axis formatting
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_yticks([0, 2, 4, 6, 8]) 
    ax.set_xlabel(f"{gene1} expression (log₂(TPM + 1))", fontsize=22, labelpad = 12)
    ax.set_ylabel(f"{gene2} expression (log₂(TPM + 1))", fontsize=22, labelpad = 12)
    ax.tick_params(labelsize=29, colors="black", width=3, length = 6)

    # Border around plot
    for spine in ax.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(0.8)

    ax.set_title(
        title_suffix,
        fontsize=27, fontweight="bold",
        loc="left",
        pad=25,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFFFCC", edgecolor="black"),
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    # Save smaller version for GitHub README
    github_path = output_path.replace(".png", "_github.png")
    plt.savefig(github_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved plot: {output_path}")


# ============================================================
# Main analysis
# ============================================================

def main():
    if not os.path.exists(EXPRESSION_FILE):
        print(f"Error: {EXPRESSION_FILE} not found. Run 01_data_prep.py first.")
        sys.exit(1)

    os.makedirs(RESULTS_DIR, exist_ok=True)

    # ==========================
    # Pan-cancer analysis
    # ==========================
    print("=== Pan-Cancer Correlation ===")

    # Load gene pair
    expr_pan = load_gene_pair(EXPRESSION_FILE, GENE1, GENE2)
    print(f"  Samples: {len(expr_pan)}")

    # Spearman correlation
    rho, pval = spearman_correlation(expr_pan[GENE1], expr_pan[GENE2])
    pval_str = format_pvalue(pval)
    print(f"  Spearman rho = {rho}, P {pval_str}")

    # TLS regression
    slope, intercept = tls_regression(expr_pan[GENE1].values, expr_pan[GENE2].values)
    print(f"  TLS: slope = {slope:.4f}, intercept = {intercept:.4f}")

    # Annotation and plot
    annotation = build_annotation(rho, pval_str)
    make_scatter_plot(
        expr_pan, GENE1, GENE2, slope, intercept, annotation,
        XLIM, YLIM,
        point_color="#FF7F50",  
        title_suffix="Pan-Cancer",
        output_path=os.path.join(RESULTS_DIR, "correlation_pancancer.png"),
    )

    # Save high-expression samples
    high_expr = expr_pan[(expr_pan[GENE1] > 6) & (expr_pan[GENE2] > 6)]
    high_expr.to_csv(os.path.join(RESULTS_DIR, f"{GENE1}_{GENE2}_top.csv"))
    print(f"  High expression samples (both > 6): {len(high_expr)}")

    # ==========================
    # AML analysis
    # ==========================
    print("\n=== AML Correlation ===")

    # Load gene pair filtered to AML samples
    expr_aml = load_gene_pair(EXPRESSION_FILE, GENE1, GENE2, sample_filter=AML_PREFIX)
    print(f"  AML samples: {len(expr_aml)}")

    # Spearman correlation
    rho_aml, pval_aml = spearman_correlation(expr_aml[GENE1], expr_aml[GENE2])
    pval_aml_str = format_pvalue(pval_aml)
    print(f"  Spearman rho = {rho_aml}, P {pval_aml_str}")

    # TLS regression
    slope_aml, intercept_aml = tls_regression(
        expr_aml[GENE1].values, expr_aml[GENE2].values
    )
    print(f"  TLS: slope = {slope_aml:.4f}, intercept = {intercept_aml:.4f}")

    # Annotation and plot
    annotation_aml = build_annotation(rho_aml, pval_aml_str)
    make_scatter_plot(
        expr_aml, GENE1, GENE2, slope_aml, intercept_aml, annotation_aml,
        XLIM, YLIM,
        point_color="#DB7093",  # R used #DB7093 for AML
        title_suffix="AML",
        output_path=os.path.join(RESULTS_DIR, "correlation_aml.png"),
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
