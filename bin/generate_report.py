#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import base64
from datetime import datetime

def image_to_base64(img_path):
    if not os.path.exists(img_path):
        return ""
    with open(img_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode('utf-8')

def generate_html(args):
    # Load data
    summary_df = pd.read_csv(args.summary_tsv, sep="\t")
    details_df = pd.read_csv(args.details_tsv, sep="\t")
    
    # Read flagstat
    flagstat_content = ""
    if os.path.exists(args.flagstat):
        with open(args.flagstat, 'r') as f:
            flagstat_content = f.read()

    # Convert plots to base64
    edit_plot_b64 = image_to_base64(args.edit_plot)
    cov_plot_b64 = image_to_base64(args.cov_plot)

    # HTML Template
    html_template = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>PolySNiP-nf Final Report - {args.sample}</title>
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; color: #333; max-width: 1200px; margin: 0 auto; padding: 20px; background-color: #f4f7f6; }}
            h1, h2, h3 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
            .card {{ background: white; padding: 20px; margin-bottom: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ text-align: left; padding: 12px; border-bottom: 1px solid #ddd; }}
            th {{ background-color: #3498db; color: white; }}
            tr:hover {{ background-color: #f5f5f5; }}
            .plot-container {{ text-align: center; margin: 20px 0; }}
            .plot-container img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; }}
            pre {{ background: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 5px; overflow-x: auto; }}
            .params-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; }}
            .stat-val {{ font-weight: bold; color: #e67e22; }}
        </style>
    </head>
    <body>
        <h1>PolySNiP-nf Analysis Report: {args.sample}</h1>
        <p>Report generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

        <div class="card">
            <h2>Pipeline Parameters</h2>
            <div class="params-grid">
                <div><strong>sgRNA Sequence:</strong> <span class="stat-val">{args.sgrna}</span></div>
                <div><strong>Quantification Window:</strong> <span class="stat-val">±{args.window}bp</span></div>
                <div><strong>Min Edit Freq:</strong> <span class="stat-val">{args.min_edit_freq}%</span></div>
                <div><strong>Min Variant Reads:</strong> <span class="stat-val">{args.min_variant_reads}</span></div>
            </div>
        </div>

        <div class="card">
            <h2>Alignment Statistics (Flagstat)</h2>
            <pre>{flagstat_content}</pre>
        </div>

        <div class="card">
            <h2>Coverage Depth Histogram</h2>
            <div class="plot-container">
                <img src="data:image/png;base64,{cov_plot_b64}" alt="Coverage Plot">
            </div>
        </div>

        <div class="card">
            <h2>Editing Frequencies</h2>
            <div class="plot-container">
                <img src="data:image/png;base64,{edit_plot_b64}" alt="Editing Plot">
            </div>
            <h3>Summary Table</h3>
            {summary_df.to_html(classes='table', index=False)}
        </div>

        <div class="card">
            <h2>Detailed Mutations (Filtered)</h2>
            <p>Showing mutations with frequency &ge; {args.min_edit_freq}%</p>
            {details_df.to_html(classes='table', index=False)}
        </div>

        <div class="card">
            <h2>Read QC (fastp) Summary</h2>
            <p>Detailed QC metrics are available in the <code>qc/</code> results directory via MultiQC.</p>
        </div>
    </body>
    </html>
    """
    
    with open(args.output, "w") as f:
        f.write(html_template)
    print(f"Final report saved to {args.output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate final HTML report.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--sgrna", required=True)
    parser.add_argument("--window", required=True)
    parser.add_argument("--min_edit_freq", required=True)
    parser.add_argument("--min_variant_reads", required=True)
    parser.add_argument("--flagstat", required=True)
    parser.add_argument("--edit_plot", required=True)
    parser.add_argument("--summary_tsv", required=True)
    parser.add_argument("--details_tsv", required=True)
    parser.add_argument("--cov_plot", required=True)
    parser.add_argument("--output", required=True)
    
    args = parser.parse_args()
    generate_html(args)
