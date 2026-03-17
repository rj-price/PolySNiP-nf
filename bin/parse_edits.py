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
        pos = record.seq.find(sgrna_seq)
        if pos != -1:
            cut_sites[record.id] = pos + sgrna_len - 3
            continue
        
        rc_sgrna = str(Seq(sgrna_seq).reverse_complement())
        pos = record.seq.find(rc_sgrna)
        if pos != -1:
            cut_sites[record.id] = pos + 3
            
    return cut_sites

def get_edit_status(bam_file, reference_fasta, cut_sites, window_size):
    """
    Parses BAM to count modified vs unmodified reads per homoeologue.
    Also returns detailed mutations and top allele sequences.
    """
    summary_results = []
    mutation_details = []
    allele_data = []
    
    ref_seqs = {record.id: str(record.seq) for record in SeqIO.parse(reference_fasta, "fasta")}
    sam_in = pysam.AlignmentFile(bam_file, "rb")
    
    for ref_id, cut_pos in cut_sites.items():
        if ref_id not in ref_seqs:
            continue
            
        win_start = cut_pos - window_size
        win_end = cut_pos + window_size
        
        unmodified = 0
        modified = 0
        mutation_counts = {}
        allele_sequences = {} # sequence -> count

        for read in sam_in.fetch(ref_id, win_start, win_end + 1):
            if read.is_unmapped: continue

            read_is_modified = False
            read_mutations = []
            
            pairs = read.get_aligned_pairs()
            
            # Find indices in pairs that correspond to our window
            window_pairs = [p for p in pairs if p[1] is not None and win_start <= p[1] <= win_end]
            
            if not window_pairs or window_pairs[0][1] > win_start or window_pairs[-1][1] < win_end:
                continue

            # Construct read's version of the window sequence
            full_window_seq = ""
            start_idx = -1
            for idx, p in enumerate(pairs):
                if p[1] == win_start:
                    start_idx = idx
                    break
            
            if start_idx == -1: continue

            curr_idx = start_idx
            while curr_idx < len(pairs):
                q_pos, r_pos = pairs[curr_idx]
                if r_pos is not None and r_pos > win_end:
                    break
                
                if q_pos is not None:
                    full_window_seq += read.query_sequence[q_pos]
                else:
                    full_window_seq += "-" # Deletion
                curr_idx += 1

            allele_sequences[full_window_seq] = allele_sequences.get(full_window_seq, 0) + 1

            # Mutation detection logic
            i = 0
            while i < len(pairs):
                q_pos, r_pos = pairs[i]
                if r_pos is not None and win_start <= r_pos <= win_end:
                    if q_pos is None:
                        del_start, del_len = r_pos, 0
                        while i < len(pairs) and pairs[i][0] is None and pairs[i][1] is not None:
                            del_len += 1
                            i += 1
                        read_mutations.append((del_start, "DEL", del_len, f"del_{del_len}bp"))
                        read_is_modified = True
                        continue
                    else:
                        rb, qb = ref_seqs[ref_id][r_pos].upper(), read.query_sequence[q_pos].upper()
                        if rb != qb:
                            read_mutations.append((r_pos, "SNP", 1, f"{rb}>{qb}"))
                            read_is_modified = True
                elif r_pos is None and q_pos is not None:
                    prev_r = pairs[i-1][1] if i > 0 else None
                    next_r = pairs[i+1][1] if i < len(pairs) - 1 else None
                    if (prev_r is not None and win_start <= prev_r <= win_end) or \
                       (next_r is not None and win_start <= next_r <= win_end):
                        ins_pos, ins_seq, ins_len = (prev_r if prev_r is not None else next_r), "", 0
                        while i < len(pairs) and pairs[i][1] is None and pairs[i][0] is not None:
                            ins_seq += read.query_sequence[pairs[i][0]]
                            ins_len += 1
                            i += 1
                        read_mutations.append((ins_pos, "INS", ins_len, f"ins_{ins_seq}"))
                        read_is_modified = True
                        continue
                i += 1
            
            if read_is_modified:
                modified += 1
                for mut in set(read_mutations):
                    mutation_counts[mut] = mutation_counts.get(mut, 0) + 1
            else:
                unmodified += 1
        
        total = modified + unmodified
        summary_results.append({
            "Homoeologue": ref_id, "Total_Reads": total, "Unmodified_Reads": unmodified,
            "Modified_Reads": modified, "Efficiency": (modified / total * 100) if total > 0 else 0
        })

        for (pos, mtype, mlen, mchange), count in mutation_counts.items():
            mutation_details.append({
                "Homoeologue": ref_id, "Position": pos, "Type": mtype, "Length": mlen,
                "Change": mchange, "Count": count, "Frequency (%)": (count / total * 100) if total > 0 else 0
            })
            
        for seq, count in allele_sequences.items():
            allele_data.append({
                "Homoeologue": ref_id, "Sequence": seq, "Count": count,
                "Frequency (%)": (count / total * 100) if total > 0 else 0
            })
        
    sam_in.close()
    return pd.DataFrame(summary_results), pd.DataFrame(mutation_details), pd.DataFrame(allele_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse BAM for CRISPR edits.")
    parser.add_argument("--vcf", help="Input VCF file (optional)")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--sgrna", required=True, help="sgRNA sequence")
    parser.add_argument("--window", type=int, default=10, help="Quantification window size")
    parser.add_argument("--min_freq", type=float, default=1.0, help="Min freq (%) to report")
    parser.add_argument("--output", required=True, help="Output summary TSV")
    parser.add_argument("--output_details", required=True, help="Output details TSV")
    parser.add_argument("--output_alleles", required=True, help="Output alleles TSV")
    
    args = parser.parse_args()
    cut_sites = find_sgrna_cut_site(args.ref, args.sgrna)
    
    if not cut_sites:
        df_summary = pd.DataFrame(columns=["Homoeologue", "Total_Reads", "Unmodified_Reads", "Modified_Reads", "Efficiency"])
        df_details = pd.DataFrame(columns=["Homoeologue", "Position", "Type", "Length", "Change", "Count", "Frequency (%)"])
        df_alleles = pd.DataFrame(columns=["Homoeologue", "Sequence", "Count", "Frequency (%)"])
    else:
        df_summary, df_details, df_alleles = get_edit_status(args.bam, args.ref, cut_sites, args.window)
        
    df_summary.to_csv(args.output, sep="\t", index=False)
    
    if not df_details.empty:
        df_details = df_details[df_details["Frequency (%)"] >= args.min_freq]
        df_details.sort_values(["Homoeologue", "Position", "Count"], ascending=[True, True, False], inplace=True)
    df_details.to_csv(args.output_details, sep="\t", index=False)

    if not df_alleles.empty:
        df_alleles = df_alleles[df_alleles["Frequency (%)"] >= args.min_freq]
        df_alleles.sort_values(["Homoeologue", "Frequency (%)"], ascending=[True, False], inplace=True)
    df_alleles.to_csv(args.output_alleles, sep="\t", index=False)
