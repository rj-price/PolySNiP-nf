#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


def plot_summaries(tsv_file, output_pdf):
    """
    Generates CRISPResso-style plots from the parsed TSV.
    """
    df = pd.read_csv(tsv_file, sep="\t")

    if df.empty:
        print("Warning: Input TSV is empty. No plots generated.")
        return

    # Set style
    sns.set_theme(style="whitegrid")

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 20))

    # Plot 1: Absolute counts (Stacked Bar Chart)
    df_plot = df.melt(
        id_vars="Homoeologue",
        value_vars=["Unmodified_Reads", "Modified_Reads"],
        var_name="Status",
        value_name="Count",
    )

    # We need to ensure Unmodified is on bottom or top as per convention
    # pivot back to wide for stacked bar
    df_wide = df.set_index("Homoeologue")[["Unmodified_Reads", "Modified_Reads"]]
    df_wide.plot(kind="bar", stacked=True, ax=ax1, color=["#7fbf7b", "#af8dc3"])

    ax1.set_title("Absolute Read Counts per Homoeologue", fontsize=16)
    ax1.set_ylabel("Number of Reads")
    ax1.set_xlabel("Homoeologue")
    ax1.legend(["Unmodified", "Modified"])

    # Plot 2: Efficiency (Percentage)
    sns.barplot(data=df, x="Homoeologue", y="Efficiency", ax=ax2, color="#91bfdb")

    ax2.set_title("Editing Efficiency per Homoeologue", fontsize=16)
    ax2.set_ylabel("Editing Efficiency (%)")
    ax2.set_xlabel("Homoeologue")
    ax2.tick_params(axis="x", rotation=90)
    ax2.set_ylim(0, 100)

    # Add labels on top of bars
    for p in ax2.patches:
        ax2.annotate(
            f"{p.get_height():.1f}%",
            (p.get_x() + p.get_width() / 2.0, p.get_height()),
            ha="center",
            va="center",
            xytext=(0, 9),
            textcoords="offset points",
        )

    plt.tight_layout()
    plt.savefig(output_pdf)
    
    # Also save as PNG for HTML report
    output_png = output_pdf.replace(".pdf", ".png")
    plt.savefig(output_png)
    
    print(f"Plots saved to {output_pdf} and {output_png}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate CRISPR editing plots.")
    parser.add_argument(
        "--input", required=True, help="Input TSV file from parse_edits.py"
    )
    parser.add_argument("--output", required=True, help="Output PDF file")

    args = parser.parse_args()

    plot_summaries(args.input, args.output)
