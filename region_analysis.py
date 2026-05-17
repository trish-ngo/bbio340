#!/usr/bin/env python3
# Comparative Sequence Conservation of Epilepsy-Associated Ion Channel Genes
# Canis lupus familiaris (Beagle) vs. Homo sapiens
# Authors: Trish Ngô 


from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
import sys
import random
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings('ignore')


# ============================================================================
# SECTION 1: SEQUENCE LOADING
# ============================================================================

def load_sequence_from_file(filename):
    """Load the first sequence from a FASTA file (.faa / .fna)."""
    try:
        records = list(SeqIO.parse(filename, "fasta"))
        if records:
            return str(records[0].seq).replace('-', '')
        print(f"  Warning: No sequences found in {filename}")
        return None
    except Exception as e:
        print(f"  Error loading sequence from {filename}: {e}")
        return None


def extract_sequence_from_blast_output(filename):
    """
    Extract the Query sequence from a BLASTP text-format output file.
    Reads lines matching 'Query  <start>  <seq_fragment>  <end>' within
    the first alignment block (stops at the second hit header).
    """
    query_re = re.compile(r'^Query\s+(\d+)\s+([A-Za-z-]+)\s+(\d+)')
    try:
        fragments = {}
        in_first_hit = False
        hit_count = 0
        with open(filename, encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                if line.startswith('>TR:') or line.startswith('>SP:'):
                    hit_count += 1
                    if hit_count == 1:
                        in_first_hit = True
                    elif hit_count == 2:
                        break
                    continue
                if not in_first_hit:
                    continue
                m = query_re.match(line)
                if m:
                    start = int(m.group(1))
                    frag = m.group(2).replace('-', '')
                    fragments[start] = frag
        if not fragments:
            print(f"  Warning: Could not extract sequence from BLAST output {filename}")
            return None
        return ''.join(f for _, f in sorted(fragments.items()))
    except Exception as e:
        print(f"  Error extracting sequence from BLAST output {filename}: {e}")
        return None


# ============================================================================
# SECTION 2: ALIGNMENT & IDENTITY
# ============================================================================

def calculate_percent_identity(seq1, seq2, sequence_type='protein'):
    """Percent identity from pairwise alignment [0, 100]."""
    if isinstance(seq1, Seq):
        seq1 = str(seq1)
    if isinstance(seq2, Seq):
        seq2 = str(seq2)
    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -2)
    if not alignments:
        return 0.0
    best = alignments[0]
    a1, a2 = best.seqA, best.seqB
    matches = 0
    total = 0
    for i in range(len(a1)):
        c1, c2 = a1[i].upper(), a2[i].upper()
        if c1 == '-' and c2 == '-':
            continue
        total += 1
        if c1 == c2 and c1 != '-':
            matches += 1
    if total == 0:
        return 0.0
    return (matches / total) * 100


def get_alignment_score(seq1, seq2, sequence_type='protein'):
    """Get the alignment score for two sequences."""
    if isinstance(seq1, Seq):
        seq1 = str(seq1)
    if isinstance(seq2, Seq):
        seq2 = str(seq2)
    alns = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -2)
    return alns[0].score if alns else 0.0


def get_aligned_part(seq1, seq2, sequence_type='protein'):
    """For long seqs, extract the part that actually aligns."""
    if isinstance(seq1, Seq):
        seq1 = str(seq1)
    if isinstance(seq2, Seq):
        seq2 = str(seq2)
    alns = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -2)
    if not alns:
        return seq1, 0.0
    aln = alns[0]
    a1, a2 = aln.seqA, aln.seqB
    part = ''.join(c1 for c1, c2 in zip(a1, a2) if c2 != '-' and c1 != '-')
    return part, aln.score


# ============================================================================
# SECTION 3: SHUFFLING STRATEGIES
# ============================================================================

def shuffle_sequence(seq):
    """Randomly permute residue order (original main.py method)."""
    s = list(str(seq))
    random.shuffle(s)
    return ''.join(s)


def shuffle_nucleotide_positions(seq, codon_aware=True):
    """
    Shuffle nucleotide positions in a DNA/CDS sequence.
    
    codon_aware=True (default): shuffle codon units (preserves codon-usage bias)
    codon_aware=False: shuffle individual nucleotide indices
    """
    seq = str(seq).upper()
    if codon_aware:
        n = (len(seq) // 3) * 3
        codons = [seq[i:i+3] for i in range(0, n, 3)]
        remainder = seq[n:]
        random.shuffle(codons)
        return ''.join(codons) + remainder
    else:
        idx = list(range(len(seq)))
        random.shuffle(idx)
        return ''.join(seq[i] for i in idx)


# ============================================================================
# SECTION 4: DOMAIN BOUNDARIES (UniProt-based)
# ============================================================================

DOMAIN_REGIONS = {
    'SCN2A': {
        'conserved': [
            (121, 399),     # Domain I   — S1-S6 + pore
            (751, 1050),    # Domain II  — S1-S6 + pore
            (1154, 1451),   # Domain III — S1-S6 + pore
            (1490, 1790),   # Domain IV  — S1-S6 + pore
        ],
        'evolving': [
            (0, 120),       # N-terminal region
            (400, 750),     # DI-DII linker
            (1051, 1153),   # DII-DIII linker
            (1452, 1489),   # DIII-DIV linker
            (1791, 2006),   # C-terminal tail
        ],
    },
    'KCNQ2': {
        'conserved': [
            (73, 323),      # S1-S6 + pore helix + selectivity filter
        ],
        'evolving': [
            (0, 72),        # N-terminal / PAS-cap domain
            (324, 660),     # C-terminal helical domain
        ],
    },
    'SCN1A': None,  # incompatible sequence types; skip region analysis
}


# ============================================================================
# SECTION 5: MONTE CARLO TESTS
# ============================================================================

def calc_pval(human_seq, dog_seq, sequence_type='protein', n=100):
    """
    Calculate p-value using Monte Carlo randomization test.
    p = (# random scores >= real) / n
    """
    if len(str(human_seq)) > 50000:
        human_part, full_sc = get_aligned_part(human_seq, dog_seq, sequence_type)
    else:
        human_part = str(human_seq)
    
    real_sc = get_alignment_score(human_part, dog_seq, sequence_type)
    
    if n <= 0:
        return {'p_value': None}
    
    count = 0
    null_scores = []
    for _ in range(n):
        shuffled = shuffle_sequence(dog_seq)
        sc = get_alignment_score(human_part, shuffled, sequence_type)
        null_scores.append(sc)
        if sc >= real_sc:
            count += 1
    
    pval = count / n
    return {
        'p_value': pval,
        'real_score': real_sc,
        'null_scores': np.array(null_scores),
        'significant': pval < 0.01,
    }


def calc_pval_position_shuffle(human_seq, dog_seq, sequence_type='protein', n=100):
    """
    Calculate p-value using nucleotide-position shuffle (new method).
    For DNA: shuffles codon positions. For protein: shuffles residues.
    """
    if len(str(human_seq)) > 50000:
        human_part, full_sc = get_aligned_part(human_seq, dog_seq, sequence_type)
    else:
        human_part = str(human_seq)
    
    real_sc = get_alignment_score(human_part, dog_seq, sequence_type)
    
    if n <= 0:
        return {'p_value': None}
    
    count = 0
    null_scores = []
    shuffle_fn = shuffle_nucleotide_positions if sequence_type.lower() == 'dna' else shuffle_sequence
    
    for _ in range(n):
        if sequence_type.lower() == 'dna':
            shuffled = shuffle_fn(dog_seq, codon_aware=True)
        else:
            shuffled = shuffle_fn(dog_seq)
        sc = get_alignment_score(human_part, shuffled, sequence_type)
        null_scores.append(sc)
        if sc >= real_sc:
            count += 1
    
    pval = count / n
    return {
        'p_value': pval,
        'real_score': real_sc,
        'null_scores': np.array(null_scores),
        'significant': pval < 0.01,
    }


def extract_region(seq, spans):
    """Concatenate subsequences at (start, end) spans (0-based, end-exclusive)."""
    return ''.join(str(seq)[max(0, s):min(e, len(str(seq)))] for s, e in spans if s < len(str(seq)))


def calc_pval_region(gene_name, human_seq, dog_seq, sequence_type='protein', n=100):
    """
    Calculate p-values separately for conserved and evolving regions.
    Returns a DataFrame with results per region.
    """
    domains = DOMAIN_REGIONS.get(gene_name.upper())
    if domains is None:
        return pd.DataFrame()
    
    rows = []
    for region_name, spans in domains.items():
        h_sub = extract_region(human_seq, spans)
        d_sub = extract_region(dog_seq, spans)
        if not h_sub or not d_sub:
            continue
        
        pct_id = calculate_percent_identity(h_sub, d_sub, sequence_type)
        real_sc = get_alignment_score(h_sub, d_sub, sequence_type)
        
        count = 0
        null_scores = []
        for _ in range(n):
            shuffled = shuffle_sequence(d_sub)
            sc = get_alignment_score(h_sub, shuffled, sequence_type)
            null_scores.append(sc)
            if sc >= real_sc:
                count += 1
        
        pval = count / n
        rows.append({
            'gene': gene_name,
            'region': region_name,
            'human_aa': len(h_sub),
            'dog_aa': len(d_sub),
            'pct_identity': round(pct_id, 2),
            'raw_score': real_sc,
            'p_value': pval,
            'significant': pval < 0.01,
            'null_scores': np.array(null_scores),
        })
    
    return pd.DataFrame(rows)


def compare_epilepsy_gene(gene_name, human_seq, dog_seq, sequence_type='protein', n=100):
    """Compare human and dog sequences for a given epilepsy gene."""
    pct_id = calculate_percent_identity(human_seq, dog_seq, sequence_type)
    pval_info = calc_pval(human_seq, dog_seq, sequence_type=sequence_type, n=n)
    return {
        'gene_name': gene_name,
        'percent_identity': pct_id,
        'p_value': pval_info['p_value'],
        'real_score': pval_info['real_score'],
        'null_scores': pval_info['null_scores'],
        'significant': pval_info['significant'],
        'human_length': len(str(human_seq).replace('-', '')),
        'dog_length': len(str(dog_seq).replace('-', ''))
    }


def compare_epilepsy_gene_position_shuffle(gene_name, human_seq, dog_seq, sequence_type='protein', n=100):
    """Compare using nucleotide-position shuffle method."""
    pct_id = calculate_percent_identity(human_seq, dog_seq, sequence_type)
    pval_info = calc_pval_position_shuffle(human_seq, dog_seq, sequence_type=sequence_type, n=n)
    return {
        'gene_name': gene_name,
        'percent_identity': pct_id,
        'p_value': pval_info['p_value'],
        'real_score': pval_info['real_score'],
        'null_scores': pval_info['null_scores'],
        'significant': pval_info['significant'],
        'human_length': len(str(human_seq).replace('-', '')),
        'dog_length': len(str(dog_seq).replace('-', ''))
    }


# ============================================================================
# SECTION 6: PLOTTING FUNCTIONS
# ============================================================================

def plot_domain_map(domain_regions, sequences_loaded):
    """Figure 0: Domain architecture map."""
    genes = [g for g, v in domain_regions.items()
             if v is not None and g in sequences_loaded]
    if not genes:
        print('No protein genes with domain maps to plot.')
        return

    colour = {'conserved': '#2a9d8f', 'evolving': '#e76f51'}
    height = {'conserved': 0.55, 'evolving': 0.35}
    seen_labels = set()

    fig, axes = plt.subplots(len(genes), 1, figsize=(13, 2.5 * len(genes)))
    if len(genes) == 1:
        axes = [axes]

    for ax, gene in zip(axes, genes):
        length = len(sequences_loaded[gene]['dog_seq'])
        ax.barh(0, length, height=0.2, color='#cccccc', zorder=1)

        for region_type, spans in domain_regions[gene].items():
            for i, (start, end) in enumerate(spans):
                lbl = region_type.capitalize() if region_type not in seen_labels else None
                ax.barh(0, end - start, left=start, height=height[region_type],
                        color=colour[region_type], alpha=0.88, label=lbl, zorder=2)
                seen_labels.add(region_type)
                mid = (start + end) / 2
                short = 'TM' if region_type == 'conserved' else 'Lnk'
                ax.text(mid, 0.35, f'{short}{i+1}', ha='center', va='bottom',
                        fontsize=7.5, color=colour[region_type], fontweight='bold')

        ax.set_xlim(0, length)
        ax.set_ylim(-0.55, 0.65)
        ax.set_yticks([])
        ax.set_xlabel('Residue position (aa)', fontsize=10)
        ax.set_title(gene, fontweight='bold', fontsize=12)
        ax.spines[['top', 'right', 'left']].set_visible(False)

    handles = [
        mpatches.Patch(color='#2a9d8f', label='Conserved: TM helices + pore'),
        mpatches.Patch(color='#e76f51', label='Evolving: linkers + terminals'),
        mpatches.Patch(color='#cccccc', label='Backbone'),
    ]
    fig.legend(handles=handles, loc='upper center', ncol=3, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.suptitle('Figure 0 — Domain Architecture\nConserved vs. Evolving Regions',
                 fontsize=11, y=1.08)
    plt.tight_layout()
    plt.savefig('fig0_domain_map.png', dpi=150, bbox_inches='tight')
    print('Saved: fig0_domain_map.png')
    plt.show()


def plot_mc_distributions(mc_dict, title_suffix='', save_path=None):
    """Figure 1: Monte Carlo null distributions."""
    n = len(mc_dict)
    fig, axes = plt.subplots(1, n, figsize=(5.5*n, 4.2), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, (gene, res) in zip(axes, mc_dict.items()):
        null = res['null_scores']
        ax.hist(null, bins=35, color='#457b9d', alpha=0.75, edgecolor='white',
                linewidth=0.4, label='Null (shuffled)')
        ax.axvline(res['real_score'], color='#e63946', linewidth=2.2, zorder=5,
                   label=f"Observed (p={res['p_value']:.3f})")
        tail = null[null >= res['real_score']]
        if len(tail):
            ax.hist(tail, bins=35, color='#e63946', alpha=0.3,
                    range=(null.min(), null.max()))
        sig = 'p<0.01 ✓' if res['significant'] else 'n.s.'
        ax.set_title(f"{gene}  [{sig}]", fontweight='bold', fontsize=12)
        ax.set_xlabel('Alignment score', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.legend(fontsize=8.5)

    fig.suptitle(f'Monte Carlo Null Distributions vs. Observed Score\n({title_suffix})',
                 fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'Saved: {save_path}')
    plt.show()


def plot_region_bars(region_df, save_path=None):
    """Figure 2: Region-aware conservation bar chart."""
    if region_df.empty:
        print('No region results to plot.')
        return

    palette = {'conserved': '#2a9d8f', 'evolving': '#e76f51'}
    labels = [f"{r['gene']}\n({r['region']})" for _, r in region_df.iterrows()]
    colours = [palette.get(r, '#aaa') for r in region_df['region']]
    x = np.arange(len(region_df))

    fig, ax = plt.subplots(figsize=(max(7, 2.8*len(region_df)), 4.8))
    ax.bar(x, region_df['pct_identity'], color=colours,
           edgecolor='white', linewidth=0.6, width=0.55)

    for i, (_, row) in enumerate(region_df.iterrows()):
        star = '***' if row['significant'] else 'n.s.'
        ax.text(i, row['pct_identity'] + 0.8, star,
                ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('% Sequence Identity (human vs. beagle)', fontsize=10)
    ax.set_ylim(0, min(108, region_df['pct_identity'].max() + 10))
    ax.set_title(
        'Figure 2 — Region-Aware Conservation: TM vs. Linker/Terminal Regions\n'
        '(*** = p<0.01 by Monte Carlo test)',
        fontsize=11, fontweight='bold')

    legend_handles = [
        mpatches.Patch(color='#2a9d8f', label='Conserved (TM helices / pore)'),
        mpatches.Patch(color='#e76f51', label='Evolving (linkers / terminal tails)'),
    ]
    ax.legend(handles=legend_handles, fontsize=9, frameon=False)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'Saved: {save_path}')
    plt.show()


# ============================================================================
# SECTION 7: MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Epilepsy Gene Sequence Conservation Calculator")
    print("Comparing Human vs. Dog (Beagle) Sequences")
    print("=" * 70)
    print()
    
    comparisons = [
        {
            'gene': 'KCNQ2',
            'human_file': 'Human_KCNQ2_datasets/ncbi_dataset/data/protein.faa',
            'dog_file': 'KCNQ2_beagle.txt',
            'sequence_type': 'protein'
        },
        {
            'gene': 'SCN1A',
            'human_file': 'Human_SCN1A_datasets/ncbi_dataset/data/gene.fna',
            'dog_file': 'SCN1A_beagle.txt',
            'sequence_type': 'dna'
        },
        {
            'gene': 'SCN2A',
            'human_file': 'Human_SCN2A_datasets/ncbi_dataset/data/protein.faa',
            'dog_file': 'SCN2A_beagle.txt',
            'sequence_type': 'protein'
        }
    ]
    
    sequences = {}
    
    for comp in comparisons:
        print(f"Loading sequences for {comp['gene']}...")
        
        human_seq = load_sequence_from_file(comp['human_file'])
        if not human_seq:
            print(f"  Warning: Could not load human sequence from {comp['human_file']}")
            continue
        
        dog_seq = extract_sequence_from_blast_output(comp['dog_file'])
        if not dog_seq:
            print(f"  Warning: Could not extract dog sequence from {comp['dog_file']}")
            continue
        
        print(f"  Human sequence length: {len(human_seq)}")
        print(f"  Dog sequence length: {len(dog_seq)}")
        
        sequences[comp['gene']] = {
            'human_seq': human_seq,
            'dog_seq': dog_seq,
            'sequence_type': comp['sequence_type']
        }
        print()
    
    if not sequences:
        print("[ERROR] No sequences loaded. Check file paths and try again.")
        return
    
    # Whole-sequence summary
    print("=" * 70)
    print("Whole-Sequence Summary")
    print("=" * 70)
    print(f"{'Gene':<10} {'% Identity':>12} {'Raw Score':>12} {'Human Len':>12} {'Dog Len':>12}")
    print("-" * 70)
    
    results = []
    for gene, seqs in sequences.items():
        pct_id = calculate_percent_identity(seqs['human_seq'], seqs['dog_seq'],
                                           seqs['sequence_type'])
        raw_score = get_alignment_score(seqs['human_seq'], seqs['dog_seq'],
                                       seqs['sequence_type'])
        results.append({
            'gene': gene,
            'pct_identity': pct_id,
            'raw_score': raw_score,
            'human_length': len(seqs['human_seq']),
            'dog_length': len(seqs['dog_seq']),
        })
        print(f"{gene:<10} {pct_id:>10.2f}%  {raw_score:>12.1f}  "
              f"{len(seqs['human_seq']):>12}  {len(seqs['dog_seq']):>12}")
    
    print()
    
    # Monte Carlo tests
    N_ITER = 500
    print("=" * 70)
    print(f"Monte Carlo Randomization Tests (n={N_ITER} iterations)")
    print("=" * 70)
    print()
    
    residue_mc = {}
    print("Residue-shuffle method (original):")
    for gene, seqs in sequences.items():
        print(f"  {gene}...", end=' ', flush=True)
        res = compare_epilepsy_gene(gene, seqs['human_seq'], seqs['dog_seq'],
                                    seqs['sequence_type'], n=N_ITER)
        residue_mc[gene] = res
        sig = '*** p<0.01' if res['significant'] else 'n.s.'
        print(f"p={res['p_value']:.4f}  {sig}")
    print()
    
    position_mc = {}
    print("Position-shuffle method (new):")
    for gene, seqs in sequences.items():
        if seqs['sequence_type'] == 'dna':
            label = 'codon-position'
        else:
            label = 'residue (protein)'
        print(f"  {gene} ({label})...", end=' ', flush=True)
        res = compare_epilepsy_gene_position_shuffle(gene, seqs['human_seq'], seqs['dog_seq'],
                                                     seqs['sequence_type'], n=N_ITER)
        position_mc[gene] = res
        sig = '*** p<0.01' if res['significant'] else 'n.s.'
        print(f"p={res['p_value']:.4f}  {sig}")
    print()
    
    # Region-aware tests
    print("=" * 70)
    print("Region-Aware Monte Carlo Tests")
    print("=" * 70)
    print()
    region_frames = []
    for gene, seqs in sequences.items():
        if seqs['sequence_type'] != 'protein':
            print(f"  [skip] {gene}: seq_type={seqs['sequence_type']} "
                  f"(domain boundaries are protein-based)")
            continue
        print(f"  {gene}...", flush=True)
        df = calc_pval_region(gene, seqs['human_seq'], seqs['dog_seq'],
                             seqs['sequence_type'], n=N_ITER)
        if not df.empty:
            region_frames.append(df)
            print(df[['gene', 'region', 'human_aa', 'dog_aa',
                     'pct_identity', 'p_value', 'significant']].to_string(index=False))
            print()
    
    region_df = pd.concat(region_frames, ignore_index=True) if region_frames else pd.DataFrame()
    print()
    
    # Figures
    print("=" * 70)
    print("Generating Figures")
    print("=" * 70)
    print()
    
    print("Figure 0: Domain architecture map...")
    plot_domain_map(DOMAIN_REGIONS, sequences)
    print()
    
    print("Figure 1a: MC null distributions (residue shuffle)...")
    plot_mc_distributions({k: {'real_score': v['real_score'], 'null_scores': v['null_scores'],
                                'p_value': v['p_value'], 'significant': v['significant']}
                          for k, v in residue_mc.items()},
                         title_suffix='residue shuffle — original method',
                         save_path='fig1a_mc_residue.png')
    print()
    
    print("Figure 1b: MC null distributions (position shuffle)...")
    plot_mc_distributions({k: {'real_score': v['real_score'], 'null_scores': v['null_scores'],
                                'p_value': v['p_value'], 'significant': v['significant']}
                          for k, v in position_mc.items()},
                         title_suffix='codon-position shuffle (DNA) / residue shuffle (protein)',
                         save_path='fig1b_mc_position.png')
    print()
    
    print("Figure 2: Region-aware conservation bar chart...")
    plot_region_bars(region_df, save_path='fig2_region_bars.png')
    print()
    
    # Final summary
    print("=" * 70)
    print("Final Summary Table")
    print("=" * 70)
    print(f"{'Gene':<10} {'% Identity':>12} {'p (residue)':>14} {'Sig (res)':>12} "
          f"{'p (position)':>14} {'Sig (pos)':>12}")
    print("-" * 70)
    
    for gene in sequences.keys():
        r = residue_mc.get(gene, {})
        pos = position_mc.get(gene, {})
        sig_res = '***' if r.get('significant', False) else 'n.s.'
        sig_pos = '***' if pos.get('significant', False) else 'n.s.'
        
        pct_id = calculate_percent_identity(sequences[gene]['human_seq'],
                                           sequences[gene]['dog_seq'],
                                           sequences[gene]['sequence_type'])
        print(f"{gene:<10} {pct_id:>10.2f}%  "
              f"{r.get('p_value', float('nan')):>13.4f}  {sig_res:>12}  "
              f"{pos.get('p_value', float('nan')):>13.4f}  {sig_pos:>12}")
    print()
    
    if not region_df.empty:
        print("Region-Aware Results:")
        print(region_df[['gene', 'region', 'pct_identity', 'p_value', 'significant']].to_string(index=False))
        print()
    
    print("=" * 70)
    print("DONE! Figures saved as PNG files in current directory.")
    print("=" * 70)


if __name__ == "__main__":
    main()
