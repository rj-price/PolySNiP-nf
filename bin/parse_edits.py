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

def get_edit_status(vcf_file, bam_file, cut_sites, window_size):
    """
    Parses VCF and BAM to count modified vs unmodified reads per homoeologue.
    """
    results = []
    
    # Load VCF to find variants in the window
    vcf_in = pysam.VariantFile(vcf_file)
    
    # For each homoeologue
    for ref_id, cut_pos in cut_sites.items():
        win_start = cut_pos - window_size
        win_end = cut_pos + window_size
        
        # Initialize counts
        unmodified = 0
        modified = 0
        
        # Use BAM to iterate through reads mapping to this homoeologue
        sam_in = pysam.AlignmentFile(bam_file, "rb")
        
        for read in sam_in.fetch(ref_id, win_start, win_end):
            # Check if this read has a variant within the window
            # We can check the CIGAR or look for VCF overlaps
            # Simplest for CRISPR amplicons: does the read have any INDEL/SNP in the window?
            
            # Use VCF to see if there are any variants for this read? 
            # Better: use the CIGAR and MD tags or just see if the read matches reference exactly in the window.
            
            # Check for any mismatches or indels in the specific window
            is_modified = False
            # Get aligned pairs for this read (read_pos, ref_pos)
            aligned_pairs = read.get_aligned_pairs(with_seq=True)
            
            for read_pos, ref_pos, ref_base in aligned_pairs:
                if ref_pos is not None and win_start <= ref_pos <= win_end:
                    # If ref_pos is in window
                    if read_pos is None: # Deletion
                        is_modified = True
                        break
                    if ref_base is not None and ref_base.islower(): # Mismatch
                        is_modified = True
                        break
                elif ref_pos is None and read_pos is not None:
                    # Insertion - we need to check if it's near the window
                    # This is trickier. Let's see if the previous/next ref_pos is in window.
                    pass

            # Count insertions by checking if they are between coordinates in the window
            if not is_modified:
                for cig_op, cig_len in read.cigartuples:
                    if cig_op == 1: # Insertion
                        # Find the reference position of this insertion
                        # (Not easily available directly from cigartuples without walking)
                        pass
            
            # Alternative: check if any VCF record overlaps this read in the window
            # This relies on the variant caller.
            
            # Let's use the VCF records that overlap the window for this ref_id
            variants_in_win = list(vcf_in.fetch(ref_id, win_start, win_end))
            
            # If any variant exists in the window, we'll check if the read supports it.
            # But CRISPR amplicons often have complex variants.
            
            # Re-evaluating: The specification asks to "Count reads supporting the reference allele (Unmodified) 
            # versus reads supporting indels/SNPs within the window (Modified)."
            
            # Let's check each read for any non-M (matches) in the window area
            # A read is modified if its alignment in the window has I, D, or X.
            
            curr_ref_pos = read.reference_start
            read_is_modified = False
            for op, length in read.cigartuples:
                # 0: M, 1: I, 2: D, 3: N, 4: S, 5: H, 6: P, 7: =, 8: X
                if op in [1, 2, 8]: # I, D, X
                    # Calculate if this occurs in the window
                    op_start = curr_ref_pos
                    op_end = curr_ref_pos + (length if op != 1 else 0)
                    
                    if not (op_end < win_start or op_start > win_end):
                        read_is_modified = True
                        break
                
                if op in [0, 2, 3, 7, 8]: # Consumes reference
                    curr_ref_pos += length
            
            if read_is_modified:
                modified += 1
            else:
                unmodified += 1
        
        sam_in.close()
        results.append({
            "Homoeologue": ref_id,
            "Total_Reads": modified + unmodified,
            "Unmodified_Reads": unmodified,
            "Modified_Reads": modified,
            "Efficiency": (modified / (modified + unmodified) * 100) if (modified + unmodified) > 0 else 0
        })
        
    return pd.DataFrame(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse VCF and BAM for CRISPR edits.")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--sgrna", required=True, help="sgRNA sequence")
    parser.add_argument("--window", type=int, default=10, help="Quantification window size")
    parser.add_argument("--output", required=True, help="Output TSV file")
    
    args = parser.parse_args()
    
    # Needs Biopython for find_sgrna_cut_site
    try:
        from Bio import SeqIO
    except ImportError:
        print("Biopython is required. Please add it to the environment.")
        sys.exit(1)
        
    cut_sites = find_sgrna_cut_site(args.ref, args.sgrna)
    
    if not cut_sites:
        print(f"Warning: sgRNA {args.sgrna} not found in reference.")
        # Create empty output
        df = pd.DataFrame(columns=["Homoeologue", "Total_Reads", "Unmodified_Reads", "Modified_Reads", "Efficiency"])
    else:
        df = get_edit_status(args.vcf, args.bam, cut_sites, args.window)
        
    df.to_csv(args.output, sep="\t", index=False)
