#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_coverage(depth_file, output_png, sample_id):
    df = pd.read_csv(depth_file, sep="\t", names=["Reference", "Position", "Depth"])
    
    if df.empty:
        print(f"Warning: {depth_file} is empty.")
        return

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 6))
    
    # Plot coverage for each reference (homoeologue)
    sns.lineplot(data=df, x="Position", y="Depth", hue="Reference")
    
    plt.title(f"Coverage Depth: {sample_id}", fontsize=16)
    plt.ylabel("Depth")
    plt.xlabel("Position in Reference")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(output_png)
    print(f"Coverage plot saved to {output_png}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot coverage depth from samtools depth.")
    parser.add_argument("--input", required=True, help="Input depth file")
    parser.add_argument("--output", required=True, help="Output PNG file")
    parser.add_argument("--sample", required=True, help="Sample ID")
    
    args = parser.parse_args()
    plot_coverage(args.input, args.output, args.sample)
