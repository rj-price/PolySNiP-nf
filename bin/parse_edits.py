#!/usr/bin/env python3
import argparse
import pysam
import pandas as pd
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def find_sgrna_cut_site(reference_fasta, sgrna_seq):
    """
    Finds the cut site (3bp upstream of PAM) for each homoeologue.
    Returns a dictionary of {ref_id: cut_site_coord}.
    """
    cut_sites = {}
    sgrna_len = len(sgrna_seq)
    
    for record in SeqIO.parse(reference_fasta, "fasta"):
        # Check forward strand
        pos = record.seq.find(sgrna_seq)
        if pos != -1:
            # Cut site is typically 3bp upstream of the end of the sgRNA (if PAM is 3')
            # Assuming sgrna_seq includes PAM or is just the guide.
            # Specification says "typically 3bp upstream of the PAM".
            # Let's assume sgrna_seq is the guide and find it. 
            # If the user provides the guide, the PAM is the next 3bp.
            # We'll calculate the cut site as the end of the guide - 3.
            # This is a bit simplistic, but we'll use the end of the match as the baseline.
            cut_sites[record.id] = pos + sgrna_len - 3
            continue
        
        # Check reverse complement
        rc_sgrna = str(Seq(sgrna_seq).reverse_complement())
        pos = record.seq.find(rc_sgrna)
        if pos != -1:
            # On reverse strand, cut site is 3bp from the start of the match (if PAM is 5')
            cut_sites[record.id] = pos + 3
            
    return cut_sites

def get_edit_status(bam_file, reference_fasta, cut_sites, window_size):
    """
    Parses BAM to count modified vs unmodified reads per homoeologue.
    Also returns a detailed list of mutations found.
    """
    summary_results = []
    mutation_details = []
    
    # Load reference sequences into memory
    ref_seqs = {record.id: str(record.seq) for record in SeqIO.parse(reference_fasta, "fasta")}
    
    # Open BAM file
    sam_in = pysam.AlignmentFile(bam_file, "rb")
    
    # For each homoeologue
    for ref_id, cut_pos in cut_sites.items():
        if ref_id not in ref_seqs:
            print(f"Warning: {ref_id} not found in reference FASTA. Skipping.")
            continue
            
        win_start = cut_pos - window_size
        win_end = cut_pos + window_size
        
        # Initialize counts
        unmodified = 0
        modified = 0
        
        # To store mutations for this homoeologue: (pos, type, length, change) -> count
        mutation_counts = {}

        # Iterate through reads mapping to this homoeologue
        for read in sam_in.fetch(ref_id, win_start, win_end + 1):
            if read.is_unmapped:
                continue

            read_is_modified = False
            read_mutations = [] # Track mutations in this specific read to avoid double counting for summary
            
            pairs = read.get_aligned_pairs()
            
            # Walk through pairs to find mutations
            i = 0
            while i < len(pairs):
                q_pos, r_pos = pairs[i]
                
                # 1. Check for SNP or Deletion
                if r_pos is not None and win_start <= r_pos <= win_end:
                    # Deletion
                    if q_pos is None:
                        del_start = r_pos
                        del_len = 0
                        while i < len(pairs) and pairs[i][0] is None and pairs[i][1] is not None:
                            del_len += 1
                            i += 1
                        
                        read_mutations.append((del_start, "DEL", del_len, f"del_{del_len}bp"))
                        read_is_modified = True
                        continue # i already advanced
                    
                    # SNP
                    else:
                        ref_base = ref_seqs[ref_id][r_pos].upper()
                        query_base = read.query_sequence[q_pos].upper()
                        if ref_base != query_base:
                            read_mutations.append((r_pos, "SNP", 1, f"{ref_base}>{query_base}"))
                            read_is_modified = True
                
                # 2. Check for Insertion
                elif r_pos is None and q_pos is not None:
                    # Find if this insertion is at/in the window
                    prev_r = pairs[i-1][1] if i > 0 else None
                    next_r = pairs[i+1][1] if i < len(pairs) - 1 else None
                    
                    if (prev_r is not None and win_start <= prev_r <= win_end) or \
                       (next_r is not None and win_start <= next_r <= win_end):
                        
                        ins_pos = prev_r if prev_r is not None else next_r
                        ins_seq = ""
                        ins_len = 0
                        while i < len(pairs) and pairs[i][1] is None and pairs[i][0] is not None:
                            ins_seq += read.query_sequence[pairs[i][0]]
                            ins_len += 1
                            i += 1
                        
                        read_mutations.append((ins_pos, "INS", ins_len, f"ins_{ins_seq}"))
                        read_is_modified = True
                        continue # i already advanced

                i += 1
            
            if read_is_modified:
                modified += 1
                # Record the unique mutations for this read in the global tally
                # (A read can have multiple mutations)
                for mut in read_mutations:
                    mutation_counts[mut] = mutation_counts.get(mut, 0) + 1
            else:
                unmodified += 1
        
        total = modified + unmodified
        summary_results.append({
            "Homoeologue": ref_id,
            "Total_Reads": total,
            "Unmodified_Reads": unmodified,
            "Modified_Reads": modified,
            "Efficiency": (modified / total * 100) if total > 0 else 0
        })

        # Add mutation details to the list
        for (pos, mtype, mlen, mchange), count in mutation_counts.items():
            mutation_details.append({
                "Homoeologue": ref_id,
                "Position": pos,
                "Type": mtype,
                "Length": mlen,
                "Change": mchange,
                "Count": count,
                "Frequency (%)": (count / total * 100) if total > 0 else 0
            })
        
    sam_in.close()
    return pd.DataFrame(summary_results), pd.DataFrame(mutation_details)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse BAM for CRISPR edits.")
    parser.add_argument("--vcf", help="Input VCF file (optional, not currently used)")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--sgrna", required=True, help="sgRNA sequence")
    parser.add_argument("--window", type=int, default=10, help="Quantification window size")
    parser.add_argument("--min_freq", type=float, default=0.5, help="Minimum frequency (%) to report a mutation in the details table")
    parser.add_argument("--output", required=True, help="Output summary TSV file")
    parser.add_argument("--output_details", required=True, help="Output detailed mutations TSV file")
    
    args = parser.parse_args()
    
    cut_sites = find_sgrna_cut_site(args.ref, args.sgrna)
    
    if not cut_sites:
        print(f"Error: sgRNA {args.sgrna} not found in reference.")
        df_summary = pd.DataFrame(columns=["Homoeologue", "Total_Reads", "Unmodified_Reads", "Modified_Reads", "Efficiency"])
        df_details = pd.DataFrame(columns=["Homoeologue", "Position", "Type", "Length", "Change", "Count", "Frequency (%)"])
    else:
        df_summary, df_details = get_edit_status(args.bam, args.ref, cut_sites, args.window)
        
    df_summary.to_csv(args.output, sep="\t", index=False)
    
    # Filter details by frequency
    if not df_details.empty:
        df_details = df_details[df_details["Frequency (%)"] >= args.min_freq]
        df_details.sort_values(["Homoeologue", "Position", "Count"], ascending=[True, True, False], inplace=True)
    
    df_details.to_csv(args.output_details, sep="\t", index=False)


