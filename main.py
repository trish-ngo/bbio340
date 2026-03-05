#!/usr/bin/env python3
# epilepsy gene sequence comparison - human vs dog

from Bio import pairwise2
from Bio.Seq import Seq
import sys
import random


def calculate_percent_identity(seq1, seq2, sequence_type='protein'):
    # percent identity from pairwise alignment
    if isinstance(seq1, Seq):
        seq1 = str(seq1)
    if isinstance(seq2, Seq):
        seq2 = str(seq2)
    if sequence_type.lower() == 'dna':
        alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -2)
    else:
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
    if isinstance(seq1, Seq):
        seq1 = str(seq1)
    if isinstance(seq2, Seq):
        seq2 = str(seq2)
    if sequence_type.lower() == 'dna':
        alns = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -2)
    else:
        alns = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -2)
    return alns[0].score if alns else 0.0


def shuffle_sequence(seq):
    s = list(str(seq))
    random.shuffle(s)
    return ''.join(s)


def get_aligned_part(seq1, seq2, sequence_type='protein'):
    # for long seqs - get the part that actually aligns (so we can run rand faster)
    if isinstance(seq1, Seq):
        seq1 = str(seq1)
    if isinstance(seq2, Seq):
        seq2 = str(seq2)
    if sequence_type.lower() == 'dna':
        alns = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -2)
    else:
        alns = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -2)
    if not alns:
        return seq1, 0.0
    aln = alns[0]
    a1, a2 = aln.seqA, aln.seqB
    part = ''.join(c1 for c1, c2 in zip(a1, a2) if c2 != '-' and c1 != '-')
    return part, aln.score


def calc_pval(human_seq, dog_seq, sequence_type='protein', n=100):
    # p = (# random scores >= real) / n
    if len(str(human_seq)) > 50000:
        human_part, full_sc = get_aligned_part(human_seq, dog_seq, sequence_type)
    else:
        human_part, full_sc = human_seq, None
    real_sc = get_alignment_score(human_part, dog_seq, sequence_type)
    if n <= 0:
        return {'p_value': None}
    count = 0
    for _ in range(n):
        shuffled = shuffle_sequence(dog_seq)
        sc = get_alignment_score(human_part, shuffled, sequence_type)
        if sc >= real_sc:
            count += 1
    pval = count / n
    return {'p_value': pval}


def compare_epilepsy_gene(gene_name, human_seq, dog_seq, sequence_type='protein', n=100):
    pct_id = calculate_percent_identity(human_seq, dog_seq, sequence_type)
    pval_info = calc_pval(human_seq, dog_seq, sequence_type=sequence_type, n=n)
    return {
        'gene_name': gene_name,
        'percent_identity': pct_id,
        'p_value': pval_info['p_value'],
        'human_length': len(str(human_seq).replace('-', '')),
        'dog_length': len(str(dog_seq).replace('-', ''))
    }


def load_sequence_from_file(filename):
    try:
        from Bio import SeqIO
        records = list(SeqIO.parse(filename, "fasta"))
        if records:
            return str(records[0].seq)
        else:
            print(f"No sequences found in {filename}")
            return None
    except Exception as e:
        print(f"Error loading sequence from {filename}: {e}")
        return None


def extract_sequence_from_blast_output(filename):
    try:
        parts_list = []
        collecting = False
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('Query') and not line.startswith('Query='):
                    parts = line.split()
                    if len(parts) >= 3:
                        seq_part = parts[2]
                        if any(c.isalpha() for c in seq_part):
                            parts_list.append(seq_part)
                            collecting = True
                elif collecting and line.startswith('>'):
                    break
        
        if parts_list:
            return ''.join(parts_list)
        else:
            print(f"Could not extract sequence from BLAST output {filename}")
            return None
    except Exception as e:
        print(f"Error extracting sequence from BLAST output {filename}: {e}")
        return None


def main():
    print("=" * 60)
    print("Epilepsy Gene Sequence Conservation Calculator")
    print("Comparing Human vs. Dog Sequences")
    print("=" * 60)
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
    
    results = []
    
    for comp in comparisons:
        print(f"Loading sequences for {comp['gene']}...")
        
        human_seq = load_sequence_from_file(comp['human_file'])
        if not human_seq:
            print(f"  Warning: Could not load human sequence from {comp['human_file']}")
            continue
        
        # Load dog sequence from BLAST output file
        dog_seq = extract_sequence_from_blast_output(comp['dog_file'])
        if not dog_seq:
            print(f"  Warning: Could not extract dog sequence from {comp['dog_file']}")
            continue
        
        print(f"  Human sequence length: {len(human_seq)}")
        print(f"  Dog sequence length: {len(dog_seq)}")
        
        result = compare_epilepsy_gene(
            comp['gene'], human_seq, dog_seq,
            sequence_type=comp['sequence_type'], n=100
        )
        
        results.append(result)
        print()
    
    # Print summary results
    print("=" * 70)
    print("Summary Results")
    print("=" * 70)
    print(f"{'Gene':<10} {'% Identity':>12} {'p-value':>12} {'Human Len':>12} {'Dog Len':>12}")
    print("-" * 70)
    
    for result in results:
        pval = result['p_value']
        if pval is None:
            pval_str = "nan"
        elif pval == 0:
            pval_str = "< 0.01"
        else:
            pval_str = f"{pval:.4f}"
        print(f"{result['gene_name']:<10} "
              f"{result['percent_identity']:>10.2f}% "
              f"{pval_str:>12} "
              f"{result['human_length']:>12} "
              f"{result['dog_length']:>12}")
    
    print("=" * 70)


if __name__ == "__main__":
    main()