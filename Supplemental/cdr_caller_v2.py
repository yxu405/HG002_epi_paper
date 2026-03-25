#!/usr/bin/env python3
"""
CDR Multi-Domain Caller v2
==========================
Identifies Centromere Dip Regions (CDRs) from phased per-CpG methylation
bedgraph files derived from Oxford Nanopore long-read sequencing, and assigns
each called domain a composite reliability score (0-100).

Background
----------
CDRs are regions within alpha-satellite arrays where CpG methylation is
substantially reduced relative to the surrounding centromeric sequence. They
mark the site of active CENP-A chromatin deposition and therefore define the
functional centromere. Within a single alpha-satellite array, the vast majority
of CpG sites are highly methylated (typically 85-95%). The CDR is characterised
by a pronounced dip below this array-wide baseline.

Because CDRs can be structurally complex -- containing multiple sub-domains
separated by partial methylation recoveries -- this caller uses a sliding-window
approach to detect all contiguous regions below a user-defined threshold,
treating each as a separate domain. A composite reliability score is then
applied to distinguish well-supported domains from marginal calls that sit
barely below the detection threshold.

Algorithm overview
------------------
1. Load per-CpG methylation values from a bedgraph file.
2. Compute the array-wide mean methylation across all CpGs in the input region.
3. Set the CDR threshold = array mean - threshold_offset (default: 10 percentage
   points below the mean).
4. Apply a sliding window (default: 10 kb window, 1 kb step) and compute the
   mean methylation within each window.
5. Identify all contiguous windows whose mean falls below the CDR threshold;
   each contiguous run becomes one candidate domain.
6. For each candidate domain, compute per-CpG statistics (mean, median, min,
   max, std, nCpG) and a composite reliability score.
7. Domains are classified HIGH_CONFIDENCE or LOW_CONFIDENCE based on whether
   their score exceeds a user-defined cutoff (default: 50).
8. The CDR envelope is reported over HIGH_CONFIDENCE domains only; the
   all-domain envelope (including LOW_CONFIDENCE) is also reported for
   reference.

Reliability score formula
-------------------------
Score = 0.50 * depth_score
      + 0.30 * cpg_score
      + 0.20 * span_score

  depth_score: captures how far the domain mean methylation falls below the
    detection threshold, normalised to the threshold offset. Domains barely
    at the threshold score near 0; domains at 1.5x the offset score 1.0.
    Formula: min((threshold - domain_mean) / threshold_offset, 1.5) / 1.5

  cpg_score: captures statistical reliability based on the number of CpG
    sites within the domain. Score saturates at 1.0 for >= 100 CpGs.
    Formula: min(nCpG / 100, 1.0)

  span_score: captures domain size as a proxy for biological relevance.
    Score saturates at 1.0 for domains >= 20 kb.
    Formula: min(span_kb / 20, 1.0)

The three components are weighted to emphasise depth (the most diagnostic
feature of a true CDR) while also rewarding domains with strong CpG support
and meaningful span. All scores are multiplied by 100 to give a 0-100 scale.

Usage
-----
    python3 cdr_caller_v2.py <input.bedgraph> [options]

Options
-------
    --window      Sliding window size in bp        (default: 10000)
    --step        Step size in bp                  (default: 1000)
    --min-cpg     Min CpGs required per window     (default: 3)
    --threshold   Percentage points below mean     (default: 10)
    --min-score   Reliability score cutoff (0-100) (default: 50)
    --outdir      Output directory
    --prefix      Output file prefix

Outputs
-------
    <prefix>_CDR_domains.bed  -- BED file with one row per domain, including
                                  reliability score and confidence label
    <prefix>_CDR_stats.txt    -- Full statistics report
    <prefix>_CDR_plot.png     -- Visualisation with HIGH/LOW domains colour-coded

Dependencies
------------
    Python >= 3.7, numpy, matplotlib
"""

import argparse
import os
import sys
from itertools import groupby
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for saving to file
import matplotlib.pyplot as plt


# =============================================================================
# 1. ARGUMENT PARSING
# =============================================================================

def parse_args():
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description='CDR Multi-Domain Caller v2 with reliability scoring.'
    )
    p.add_argument('bedgraph',
                   help='Input bedgraph file: chrom, start, end, methylation%%')
    p.add_argument('--window', type=int, default=10000,
                   help='Sliding window size in bp (default: 10000)')
    p.add_argument('--step', type=int, default=1000,
                   help='Step size between windows in bp (default: 1000)')
    p.add_argument('--min-cpg', type=int, default=3,
                   help='Min CpGs required in a window to compute its mean (default: 3)')
    p.add_argument('--threshold', type=float, default=10.0,
                   help='Threshold offset: CDR threshold = array_mean - this value (default: 10.0)')
    p.add_argument('--min-score', type=float, default=50.0,
                   help='Reliability score cutoff 0-100; domains below are LOW_CONFIDENCE (default: 50)')
    p.add_argument('--outdir', type=str, default=None,
                   help='Output directory (default: same as input file)')
    p.add_argument('--prefix', type=str, default=None,
                   help='Output filename prefix (default: input filename stem)')
    return p.parse_args()


# =============================================================================
# 2. DATA LOADING
# =============================================================================

def load_bedgraph(path):
    """
    Load a 4-column bedgraph file into numpy arrays.

    Expects columns: chrom, start, end, methylation_percent
    Lines beginning with 'track', '#', or empty lines are skipped.

    Parameters
    ----------
    path : str
        Path to the bedgraph file.

    Returns
    -------
    chrom : str
        Chromosome/haplotype name (taken from the first data line).
    positions : np.ndarray
        Start positions of each CpG (integers).
    values : np.ndarray
        Methylation percentages (floats, 0-100).
    """
    positions, values, chroms = [], [], []
    with open(path) as f:
        for line in f:
            # Skip header and comment lines
            if line.startswith('track') or line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chroms.append(parts[0])
            positions.append(int(parts[1]))
            values.append(float(parts[3]))
    if not positions:
        sys.exit('ERROR: No valid records found in %s' % path)
    return chroms[0], np.array(positions), np.array(values)


# =============================================================================
# 3. CORE CDR DETECTION
# =============================================================================

def sliding_window(positions, values, window_size, step, min_cpg):
    """
    Apply a sliding window across the array and compute the mean methylation
    in each window.

    Windows are placed at every `step` bp starting from the first CpG position.
    Windows containing fewer than `min_cpg` CpGs are skipped to avoid
    unreliable means in sparsely covered regions.

    Parameters
    ----------
    positions : np.ndarray
        CpG start positions.
    values : np.ndarray
        CpG methylation values.
    window_size : int
        Width of the sliding window in bp.
    step : int
        Step size between consecutive windows in bp.
    min_cpg : int
        Minimum CpGs required to compute a window mean.

    Returns
    -------
    centers : np.ndarray
        Genomic midpoints of each valid window.
    means : np.ndarray
        Mean methylation within each valid window.
    """
    # Generate start positions for each window across the array
    coords = np.arange(positions[0], positions[-1] - window_size + 1, step)
    means, centers = [], []
    for c in coords:
        # Select CpGs falling within this window [c, c + window_size)
        mask = (positions >= c) & (positions < c + window_size)
        if mask.sum() >= min_cpg:
            means.append(np.mean(values[mask]))
            centers.append(c + window_size // 2)  # midpoint of window
    return np.array(centers), np.array(means)


def call_domains(window_centers, window_means, threshold):
    """
    Identify CDR domains as contiguous runs of windows below the threshold.

    Uses Python's itertools.groupby to find runs of True values in the
    boolean mask, where True = window mean < threshold. Each contiguous
    run becomes one candidate domain, defined by the genomic positions of
    the first and last windows in the run.

    Parameters
    ----------
    window_centers : np.ndarray
        Midpoint positions of sliding windows.
    window_means : np.ndarray
        Mean methylation of each window.
    threshold : float
        CDR detection threshold (array mean - threshold_offset).

    Returns
    -------
    domains : list of (int, int)
        List of (start, end) genomic coordinate pairs for each candidate domain.
    """
    # Boolean mask: True where the window dips below the CDR threshold
    dip_mask = window_means < threshold
    domains = []
    for key, group in groupby(enumerate(dip_mask), key=lambda x: x[1]):
        items = list(group)
        if key:  # key=True means this run is below threshold
            s, e = items[0][0], items[-1][0]
            domains.append((int(window_centers[s]), int(window_centers[e])))
    return domains


def domain_stats(positions, values, s, e):
    """
    Compute per-CpG summary statistics for a single domain.

    All individual CpG methylation values within [s, e] are collected
    and summarised. Returns None if no CpGs fall within the interval
    (degenerate domain, filtered downstream).

    Parameters
    ----------
    positions : np.ndarray
        CpG positions across the full array.
    values : np.ndarray
        CpG methylation values.
    s, e : int
        Start and end coordinates of the domain.

    Returns
    -------
    dict or None
        Dictionary of summary statistics, or None if no CpGs found.
    """
    mask = (positions >= s) & (positions <= e)
    v = values[mask]
    if len(v) == 0:
        return None
    return {
        'n_cpg':  int(mask.sum()),
        'mean':   float(np.mean(v)),
        'median': float(np.median(v)),
        'min':    float(np.min(v)),
        'max':    float(np.max(v)),
        'std':    float(np.std(v)),
    }


# =============================================================================
# 4. RELIABILITY SCORING
# =============================================================================

def reliability_score(st, threshold, threshold_offset):
    """
    Compute a composite reliability score (0-100) for a CDR domain.

    The score captures three orthogonal properties of domain quality:

    1. Depth (weight 50%): How far the domain mean methylation falls below
       the CDR threshold, normalised to the threshold offset. A domain at
       exactly the threshold scores 0; one at 1.5x the offset scores 1.0.
       This is the most diagnostically important component -- shallow domains
       that barely cross the threshold are penalised heavily.

    2. CpG count (weight 30%): Number of individual CpG observations within
       the domain, capturing statistical reliability. Saturates at 100 CpGs
       (score = 1.0). Domains with very few CpGs (<20) are penalised.

    3. Span (weight 20%): Domain size in kb as a proxy for biological
       relevance. CDRs are typically tens to hundreds of kb; very small
       domains (<5 kb) are penalised. Saturates at 20 kb.

    Parameters
    ----------
    st : dict
        Domain statistics dictionary from domain_stats().
    threshold : float
        CDR detection threshold (array mean - offset).
    threshold_offset : float
        The offset used to set the threshold (e.g., 10.0 percentage points).

    Returns
    -------
    float
        Reliability score in range [0, 100].
    """
    span_kb = st['span_kb']

    # Component 1: Depth below threshold
    # Normalise depth to the threshold offset; cap at 1.5x to avoid
    # extreme outliers dominating the score
    depth = min((threshold - st['mean']) / threshold_offset, 1.5) / 1.5

    # Component 2: CpG count (statistical power)
    # Saturates at 100 CpGs
    cpg = min(st['n_cpg'] / 100.0, 1.0)

    # Component 3: Domain span (size reliability)
    # Saturates at 20 kb
    span = min(span_kb / 20.0, 1.0)

    # Weighted composite score scaled to 0-100
    score = (0.50 * depth + 0.30 * cpg + 0.20 * span) * 100.0

    # Clamp to zero (depth can be negative if domain mean exceeds threshold,
    # which can occur when a degenerate sliding-window call captures a region
    # whose raw CpG mean is actually above threshold)
    return round(max(score, 0.0), 1)


# =============================================================================
# 5. OUTPUT WRITERS
# =============================================================================

def write_bed(domains, stats_list, chrom, out_path):
    """
    Write CDR domains to a BED-format file.

    Columns: chrom, start, end, name, reliability_score, strand,
             n_cpg, mean_meth, min_meth, max_meth, confidence

    The reliability_score is written in the BED score column (col 5),
    and a HIGH_CONFIDENCE / LOW_CONFIDENCE label is appended as col 11.
    This file is compatible with UCSC Genome Browser custom tracks.

    Parameters
    ----------
    domains : list of (int, int)
        Domain coordinate pairs.
    stats_list : list of dict
        Per-domain statistics including 'score' and 'min_score'.
    chrom : str
        Chromosome/haplotype name.
    out_path : str
        Output file path.
    """
    with open(out_path, 'w') as f:
        f.write('# CDR domains BED - reliability scored\n')
        f.write('# chrom\tstart\tend\tname\treliability_score\tstrand'
                '\tn_cpg\tmean_meth\tmin_meth\tmax_meth\tconfidence\n')
        for i, (s, e) in enumerate(domains):
            st = stats_list[i]
            conf = 'HIGH_CONFIDENCE' if st['score'] >= st['min_score'] else 'LOW_CONFIDENCE'
            f.write('%s\t%d\t%d\tCDR_domain_%d\t%.1f\t.\t%d\t%.2f\t%.2f\t%.2f\t%s\n' % (
                chrom, s, e, i + 1, st['score'], st['n_cpg'],
                st['mean'], st['min'], st['max'], conf))


def write_stats(domains, stats_list, chrom, array_mean, threshold,
                window_size, step, min_cpg, threshold_offset, min_score,
                input_file, out_path):
    """
    Write a full human-readable statistics report.

    Includes run parameters, array-wide summary, envelope coordinates
    (HIGH_CONFIDENCE domains only, plus all-domain for reference), and
    a per-domain table with all statistics and reliability scores.

    The HIGH_CONFIDENCE envelope spans from the start of the first
    HIGH_CONFIDENCE domain to the end of the last, providing a biologically
    meaningful CDR extent that excludes marginal noise calls at the periphery.

    Parameters
    ----------
    (see main() for parameter descriptions)
    """
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    n_hi = sum(1 for st in stats_list if st['score'] >= min_score)
    n_lo = len(stats_list) - n_hi

    with open(out_path, 'w') as f:
        f.write('=' * 70 + '\n')
        f.write('CDR MULTI-DOMAIN CALLER v2 - RELIABILITY-SCORED REPORT\n')
        f.write('=' * 70 + '\n')
        f.write('Run date:              %s\n' % now)
        f.write('Input file:            %s\n' % input_file)

        f.write('\n--- Analysis Parameters ---\n')
        f.write('Sliding window size:   %d bp\n' % window_size)
        f.write('Step size:             %d bp\n' % step)
        f.write('Min CpGs per window:   %d\n' % min_cpg)
        f.write('CDR threshold offset:  %.1f%%\n' % threshold_offset)
        f.write('Reliability cutoff:    %.0f / 100\n' % min_score)

        f.write('\n--- Reliability Score Formula ---\n')
        f.write('  Score = 0.50 * depth_score  (depth below threshold)\n')
        f.write('        + 0.30 * cpg_score    (nCpG, saturates at 100 CpGs)\n')
        f.write('        + 0.20 * span_score   (span, saturates at 20 kb)\n')
        f.write('  Domains with Score < %.0f are flagged LOW_CONFIDENCE\n' % min_score)

        f.write('\n--- Array-Wide Summary ---\n')
        f.write('Chromosome/region:     %s\n' % chrom)
        f.write('Array mean meth:       %.2f%%\n' % array_mean)
        f.write('CDR threshold:         %.2f%%\n' % threshold)
        f.write('Total domains called:  %d\n' % len(domains))
        f.write('  HIGH_CONFIDENCE:     %d\n' % n_hi)
        f.write('  LOW_CONFIDENCE:      %d\n' % n_lo)

        # HIGH_CONFIDENCE envelope: spans only confirmed domains
        hi_domains = [(s, e) for (s, e), st in zip(domains, stats_list)
                      if st['score'] >= min_score]
        if len(hi_domains) >= 1:
            env_start = hi_domains[0][0]
            env_end   = hi_domains[-1][1]
            env_span  = env_end - env_start
            f.write('HIGH_CONF CDR envelope: %d - %d (%d bp / %.1f kb)  [%d HIGH_CONFIDENCE domain(s)]\n' % (
                env_start, env_end, env_span, env_span / 1000, len(hi_domains)))

        # All-domain envelope reported for reference when LOW_CONFIDENCE domains exist
        if len(domains) > len(hi_domains):
            f.write('All-domain envelope:    %d - %d (%d bp / %.1f kb)  [all %d domain(s), incl. LOW_CONFIDENCE]\n' % (
                domains[0][0], domains[-1][1],
                domains[-1][1] - domains[0][0],
                (domains[-1][1] - domains[0][0]) / 1000,
                len(domains)))

        # Per-domain table
        f.write('\n--- Per-Domain Statistics ---\n')
        hdr = '%-10s %-14s %-14s %-10s %-8s %-8s %-8s %-8s %-8s %-8s %-10s %-18s\n'
        f.write(hdr % ('Domain', 'Start', 'End', 'Span(kb)', 'nCpG',
                       'Mean%', 'Med%', 'Min%', 'Max%', 'Std%', 'Score', 'Confidence'))
        f.write('-' * 120 + '\n')
        row = '%-10s %-14d %-14d %-10.1f %-8d %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-10.1f %-18s\n'
        for i, ((s, e), st) in enumerate(zip(domains, stats_list)):
            conf = 'HIGH_CONFIDENCE' if st['score'] >= min_score else 'LOW_CONFIDENCE'
            f.write(row % ('Domain_%d' % (i + 1), s, e, st['span_kb'],
                           st['n_cpg'], st['mean'], st['median'],
                           st['min'], st['max'], st['std'],
                           st['score'], conf))
        f.write('\n')


def write_plot(positions, values, window_centers, window_means,
               domains, stats_list, chrom, array_mean, threshold,
               min_score, out_path):
    """
    Generate and save a visualisation of CDR domain calls.

    The plot shows:
      - Individual CpG methylation values as a scatter plot (light blue dots)
      - The 10 kb sliding window mean as a line (navy)
      - The CDR detection threshold as a dashed horizontal line (orange)
      - The array-wide mean as a dotted horizontal line (grey)
      - Each CDR domain as a vertical shaded span:
            HIGH_CONFIDENCE domains: coloured spans (distinct colour per domain)
            LOW_CONFIDENCE domains:  grey spans
      - The HIGH_CONFIDENCE CDR envelope as a double-headed arrow annotation
        spanning from the start of the first to the end of the last
        HIGH_CONFIDENCE domain.

    Parameters
    ----------
    (see main() for parameter descriptions)
    """
    # Colour palette for HIGH_CONFIDENCE domains (cycles if >7 domains)
    hi_colors = ['#d62728', '#ff7f0e', '#9467bd', '#e377c2',
                 '#8c564b', '#2ca02c', '#17becf']
    lo_color  = '#aaaaaa'  # grey for LOW_CONFIDENCE

    fig, ax = plt.subplots(figsize=(14, 5))

    # Layer 1: individual CpG methylation values as background scatter
    ax.scatter(positions / 1e6, values, s=1, alpha=0.15, color='steelblue',
               label='Individual CpGs', zorder=1)

    # Layer 2: sliding window mean as a summary trace
    ax.plot(window_centers / 1e6, window_means, color='navy', linewidth=1.2,
            label='10 kb rolling mean', zorder=2)

    # Reference lines: CDR threshold and array-wide mean
    ax.axhline(threshold, color='darkorange', linestyle='--', linewidth=1.3,
               label='CDR threshold (%.1f%%)' % threshold, zorder=3)
    ax.axhline(array_mean, color='gray', linestyle=':', linewidth=1.1,
               label='Array mean (%.1f%%)' % array_mean, zorder=3)

    # Layer 3: shade each domain; colour = confidence level
    hi_idx = 0
    for i, ((s, e), st) in enumerate(zip(domains, stats_list)):
        is_hi = st['score'] >= min_score
        if is_hi:
            color = hi_colors[hi_idx % len(hi_colors)]
            hi_idx += 1
            alpha = 0.25
            lbl = 'Domain %d - score %.0f [HIGH] (%.3f-%.3f Mb)' % (
                i + 1, st['score'], s / 1e6, e / 1e6)
        else:
            color = lo_color
            alpha = 0.30
            lbl = 'Domain %d - score %.0f [LOW] (%.3f-%.3f Mb)' % (
                i + 1, st['score'], s / 1e6, e / 1e6)
        ax.axvspan(s / 1e6, e / 1e6, color=color, alpha=alpha,
                   label=lbl, zorder=2)

    # Layer 4: HIGH_CONFIDENCE envelope bracket annotation
    # Spans from start of first to end of last HIGH_CONFIDENCE domain,
    # drawn as a double-headed arrow near the bottom of the plot
    hi_domains = [(s, e) for (s, e), st in zip(domains, stats_list)
                  if st['score'] >= min_score]
    if hi_domains:
        env_start = hi_domains[0][0] / 1e6
        env_end   = hi_domains[-1][1] / 1e6
        env_y = 4.5  # y position for the bracket (near bottom of 0-105 axis)
        ax.annotate('', xy=(env_end, env_y), xytext=(env_start, env_y),
                    arrowprops=dict(arrowstyle='<->', color='black', lw=1.5),
                    zorder=5)
        ax.text((env_start + env_end) / 2, env_y + 1.8,
                'HIGH_CONF envelope: %.3f-%.3f Mb (%.1f kb)' % (
                    env_start, env_end,
                    (hi_domains[-1][1] - hi_domains[0][0]) / 1000),
                ha='center', va='bottom', fontsize=7.5, color='black',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.7))

    ax.set_xlabel('%s position (Mb)' % chrom, fontsize=12)
    ax.set_ylabel('CpG Methylation (%)', fontsize=12)
    ax.set_title('%s - CDR Domain Calls with Reliability Scores' % chrom,
                 fontsize=13, fontweight='bold')
    ax.set_ylim(0, 105)
    ax.legend(loc='lower left', fontsize=8, markerscale=5)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches='tight')
    plt.close()


# =============================================================================
# 6. MAIN PIPELINE
# =============================================================================

def main():
    """
    Main pipeline: load data, detect CDR domains, score, and write outputs.
    """
    args = parse_args()

    # Resolve output paths
    infile = args.bedgraph
    outdir = args.outdir if args.outdir else os.path.dirname(os.path.abspath(infile))
    prefix = args.prefix if args.prefix else os.path.splitext(os.path.basename(infile))[0]
    os.makedirs(outdir, exist_ok=True)

    bed_path   = os.path.join(outdir, prefix + '_CDR_domains.bed')
    stats_path = os.path.join(outdir, prefix + '_CDR_stats.txt')
    plot_path  = os.path.join(outdir, prefix + '_CDR_plot.png')

    # Step 1: Load per-CpG methylation data
    print('[1/4] Loading: %s' % infile)
    chrom, positions, values = load_bedgraph(infile)
    print('      %d CpGs  |  %s:%d-%d' % (len(values), chrom, positions[0], positions[-1]))

    # Step 2: Compute array-wide mean and set CDR detection threshold
    array_mean = float(np.mean(values))
    threshold  = array_mean - args.threshold
    print('[2/4] Array mean: %.2f%%  |  Threshold: %.2f%%' % (array_mean, threshold))

    # Step 3: Apply sliding window to smooth the methylation signal
    window_centers, window_means = sliding_window(
        positions, values, args.window, args.step, args.min_cpg)

    # Step 4: Call candidate domains as contiguous below-threshold windows
    domains = call_domains(window_centers, window_means, threshold)
    if not domains:
        sys.exit('No CDR domains found.')

    # Step 5: Filter degenerate domains (no CpGs in the raw data, e.g. gaps)
    domains = [(s, e) for s, e in domains
               if ((positions >= s) & (positions <= e)).sum() > 0]

    # Step 6: Compute per-domain statistics and reliability scores
    stats_list = []
    for s, e in domains:
        st = domain_stats(positions, values, s, e)
        if st is None:
            continue
        st['span_kb']   = (e - s) / 1000.0
        st['score']     = reliability_score(st, threshold, args.threshold)
        st['min_score'] = args.min_score  # store cutoff for use in writers
        stats_list.append(st)

    # Keep domains and stats aligned after any None-filtering
    domains = domains[:len(stats_list)]

    # Print domain summary to console
    print('      %d domain(s) called:' % len(domains))
    for i, ((s, e), st) in enumerate(zip(domains, stats_list)):
        conf = 'HIGH' if st['score'] >= args.min_score else 'LOW '
        print('        Domain %d [%s score=%.0f]: %s:%d-%d  %.1f kb  mean=%.2f%%  nCpG=%d' % (
            i + 1, conf, st['score'], chrom, s, e,
            st['span_kb'], st['mean'], st['n_cpg']))

    # Step 7: Write outputs
    print('[3/4] Writing outputs...')
    write_bed(domains, stats_list, chrom, bed_path)
    write_stats(domains, stats_list, chrom, array_mean, threshold,
                args.window, args.step, args.min_cpg, args.threshold,
                args.min_score, os.path.basename(infile), stats_path)
    print('      BED:   %s' % bed_path)
    print('      Stats: %s' % stats_path)

    # Step 8: Generate visualisation
    print('[4/4] Generating plot...')
    write_plot(positions, values, window_centers, window_means,
               domains, stats_list, chrom, array_mean, threshold,
               args.min_score, plot_path)
    print('      Plot:  %s' % plot_path)
    print('Done.')


if __name__ == '__main__':
    main()
