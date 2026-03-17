#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

def get_alignment_grid(ref_seq, alleles):
    """
    Creates a synchronized alignment grid for a reference and multiple alleles.
    Returns (expanded_ref, aligned_alleles_list)
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -1

    # 1. Align each allele to ref and track insertions
    # insertion_map: ref_idx -> list of inserted sequences from different alleles
    insertion_map = {i: [] for i in range(len(ref_seq) + 1)}
    
    allele_alignments = []
    for allele_seq in alleles:
        # Remove any existing gaps from the sequence for a fresh alignment
        clean_allele = allele_seq.replace('-', '')
        alignments = aligner.align(ref_seq, clean_allele)
        best_aln = alignments[0]
        
        # Extract target (ref) and query (allele) aligned strings
        # Format: (aligned_ref, aligned_query)
        ref_aln_str, query_aln_str = format(best_aln).splitlines()[:2]
        # Biopython format might include match lines, we want the raw strings
        # Actually, format(best_aln) output depends on version. 
        # Safer to use aligned property:
        target_indices, query_indices = best_aln.aligned
        
        # Build alignment strings manually from indices
        t_str = ""
        q_str = ""
        last_t = 0
        last_q = 0
        for t_range, q_range in zip(target_indices, query_indices):
            # Deletion in query (gap in query)
            while last_t < t_range[0]:
                t_str += ref_seq[last_t]
                q_str += "-"
                last_t += 1
            # Insertion in query (gap in target)
            while last_q < q_range[0]:
                t_str += "-"
                q_str += clean_allele[last_q]
                last_q += 1
            # Match/Mismatch
            while last_t < t_range[1]:
                t_str += ref_seq[last_t]
                q_str += clean_allele[last_q]
                last_t += 1
                last_q += 1
        
        # Handle trailing
        while last_t < len(ref_seq):
            t_str += ref_seq[last_t]
            q_str += "-"
            last_t += 1
        while last_q < len(clean_allele):
            t_str += "-"
            q_str += clean_allele[last_q]
            last_q += 1
            
        allele_alignments.append((t_str, q_str))

    # 2. Consolidate all insertions into a master reference grid
    # We walk through all alignments and find where gaps were inserted in target
    master_ref = ""
    # This is tricky. Let's simplify: 
    # Just use the alignment that resulted in the longest string as a base? No.
    
    # Better: Use a pointer-based merge
    # For every position in ref_seq, what is the max number of '-' inserted before it?
    max_ins_at_pos = [0] * (len(ref_seq) + 1)
    for t_str, q_str in allele_alignments:
        curr_ref_idx = 0
        curr_ins_count = 0
        for char in t_str:
            if char == '-':
                curr_ins_count += 1
            else:
                max_ins_at_pos[curr_ref_idx] = max(max_ins_at_pos[curr_ref_idx], curr_ins_count)
                curr_ref_idx += 1
                curr_ins_count = 0
        max_ins_at_pos[curr_ref_idx] = max(max_ins_at_pos[curr_ref_idx], curr_ins_count)

    # 3. Build expanded ref and re-align alleles to it
    expanded_ref = ""
    for i, char in enumerate(ref_seq):
        expanded_ref += ("-" * max_ins_at_pos[i]) + char
    expanded_ref += ("-" * max_ins_at_pos[len(ref_seq)])

    # 4. Final alignment of each allele to expanded_ref
    final_alleles = []
    for allele_seq in alleles:
        clean_allele = allele_seq.replace('-', '')
        # Align to expanded_ref, but preserve expanded_ref gaps
        # We can just use the indices from previous alignment and map them to expanded_ref coordinates
        # But a simple re-alignment with high gap penalty for expanded_ref gaps is easier
        alignments = aligner.align(expanded_ref, clean_allele)
        best_aln = alignments[0]
        
        # Extract query string corresponding to expanded_ref
        target_indices, query_indices = best_aln.aligned
        q_final = ["-"] * len(expanded_ref)
        for t_range, q_range in zip(target_indices, query_indices):
            for i in range(t_range[1] - t_range[0]):
                q_final[t_range[0] + i] = clean_allele[q_range[0] + i]
        
        final_alleles.append("".join(q_final))

    return expanded_ref, final_alleles

def plot_allele_view(alleles_tsv, ref_fasta, sgrna, output_png, window_size):
    df = pd.read_csv(alleles_tsv, sep="\t")
    if df.empty:
        print("Warning: No allele data found.")
        return

    ref_seqs = {record.id: str(record.seq) for record in SeqIO.parse(ref_fasta, "fasta")}
    homoeologues = df['Homoeologue'].unique()
    num_hom = len(homoeologues)
    
    fig, axes = plt.subplots(num_hom, 1, figsize=(18, 5 * num_hom), squeeze=False)
    
    for idx, hom in enumerate(homoeologues):
        ax = axes[idx, 0]
        hom_df = df[df['Homoeologue'] == hom].head(10)
        
        # 1. Get reference window
        ref_seq = ref_seqs[hom]
        sgrna_len = len(sgrna)
        pos = ref_seq.find(sgrna)
        is_rc = False
        if pos == -1:
            rc_sgrna = str(Seq(sgrna).reverse_complement())
            pos = ref_seq.find(rc_sgrna)
            is_rc = True
        
        cut_site_in_ref = (pos + sgrna_len - 3) if not is_rc else (pos + 3)
        win_start = cut_site_in_ref - window_size
        win_end = cut_site_in_ref + window_size
        ref_win_seq = ref_seq[win_start : win_end + 1]
        
        # 2. Synchronize alignment
        top_allele_seqs = hom_df['Sequence'].tolist()
        expanded_ref, aligned_alleles = get_alignment_grid(ref_win_seq, top_allele_seqs)
        
        # 3. Plot
        y_pos = 0
        char_width = 0.10
        
        # Plot Reference
        ax.text(-0.8, y_pos, "Reference", verticalalignment='center', fontweight='bold', fontsize=11)
        
        # Calculate cut site position in expanded ref
        # It's after the 'window_size'th non-gap character
        ref_chars_seen = 0
        cut_site_x = 0
        for i, char in enumerate(expanded_ref):
            color = 'black'
            if char != '-':
                # This is a real ref base. Check if it's part of sgRNA.
                curr_ref_pos = win_start + ref_chars_seen
                if pos <= curr_ref_pos < pos + sgrna_len:
                    color = 'red'
                
                if ref_chars_seen == window_size:
                    cut_site_x = i
                ref_chars_seen += 1
            
            ax.text(i * char_width, y_pos, char, color=color, family='monospace', fontsize=11, verticalalignment='center')
        
        # Draw cut site line
        ax.axvline(x=(cut_site_x - 0.5 + (1 if is_rc else 0)) * char_width, 
                   color='blue', linestyle='--', alpha=0.4)
        
        y_pos -= 1
        
        # Plot Alleles
        for i, (aligned_seq, (_, row)) in enumerate(zip(aligned_alleles, hom_df.iterrows())):
            freq = row['Frequency (%)']
            ax.text(-0.8, y_pos, f"{freq:.1f}%", verticalalignment='center', fontsize=9)
            
            for j, char in enumerate(aligned_seq):
                ref_char = expanded_ref[j]
                
                color = 'gray' # Match
                if char == '-':
                    color = 'black' # Deletion
                elif ref_char == '-':
                    color = 'green' # Insertion
                elif char.upper() != ref_char.upper():
                    color = 'orange' # SNP
                
                ax.text(j * char_width, y_pos, char, color=color, family='monospace', fontsize=11, verticalalignment='center')
            
            y_pos -= 1

        ax.set_title(f"Allele View: {hom}", fontsize=14, pad=20)
        ax.axis('off')
        ax.set_xlim(-1.0, len(expanded_ref) * char_width + 0.5)
        ax.set_ylim(y_pos, 1)

    plt.tight_layout()
    plt.savefig(output_png, bbox_inches='tight', dpi=200)
    print(f"Improved allele view saved to {output_png}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--sgrna", required=True)
    parser.add_argument("--window", type=int, required=True)
    parser.add_argument("--output", required=True)
    
    args = parser.parse_args()
    plot_allele_view(args.input, args.ref, args.sgrna, args.output, args.window)
