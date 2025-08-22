os.chdir('/root/Cycle/HarmoCycle_v1/Codes/')
import scanpy as sc
from HarmoCycle_v4 import * 
from HarmoCycle_utils import *

# temp_res = adata1.var.head(200).copy()

ORGANISM = 'Human'
enrichment_results = run_comprehensive_enrichment(gene_list=temp_res.index.tolist(), organism=ORGANISM)

Res1_GO_BP = enrichment_results.get('GO_Biological_Process_2021', pd.DataFrame())
plot_enrichment(
        enrichment_df=Res1_GO_BP,
        title='U5. Top GO_BP Pathways (Human)',
        top_n=10, # 最多显示15个
        cutoff=0.1
    )

enrichment_results = enrichment_results['GO_Biological_Process_2021']

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import textwrap
from matplotlib.colors import LinearSegmentedColormap # Import this new module

# --- 1. Simulate and Prepare Your Data ---
# This part remains the same. In your real script, you'll load your data here.
# For example:
# enrichment_results = pd.read_csv('your_enrichment_results.csv')
# enrichment_results = enrichment_results['GO_Biological_Process_2021'] # Uncomment if needed

# Using the simulated data for this example



# --- 2. Process Data for Plotting ---
# This part is identical to the previous version.
top_n = 10
adj_p_cutoff = 0.1
significant_results = enrichment_results[enrichment_results['Adjusted P-value'] < adj_p_cutoff].copy()
significant_results['-log10(AdjP)'] = -np.log10(significant_results['Adjusted P-value'])
significant_results = significant_results.sort_values(by='-log10(AdjP)', ascending=True)
plot_data = significant_results.tail(top_n)
plot_data['Term'] = plot_data['Term'].str.replace(r'\s\(GO:\d+\)$', '', regex=True)
plot_data['Term'] = plot_data['Term'].apply(
    lambda x: textwrap.fill(x.capitalize(), width=50, subsequent_indent='  ')
)

print("--- Data Prepared for Plotting ---")
print(plot_data[['Term', '-log10(AdjP)']])
print("\n")


# --- 3. Create the "Nature Style" Plot with Gradient ---
# *** THIS SECTION CONTAINS THE KEY CHANGES ***

# Define the gradient colors based on your palette
color_start = "#E3A897"    # Lightest color (for least significant)
color_end = "#C85B43"      # Darkest color (for most significant)
text_color = "#333333"
grid_color = "#D3D3D3"

# Create a custom colormap for the gradient
cmap = LinearSegmentedColormap.from_list("custom_gradient", [color_start, color_end])

# Generate a list of colors from the colormap
# One color for each bar, creating the gradient effect
n_bars = len(plot_data)
bar_colors = cmap(np.linspace(0, 1, n_bars))

# Create the figure and axes
fig, ax = plt.subplots(figsize=(8, 6))

# Create the horizontal bars, now using the list of gradient colors
bars = ax.barh(
    plot_data['Term'],
    plot_data['-log10(AdjP)'],
    color=bar_colors, # <-- The main change is here!
    height=0.7
)

# Add text labels (Overlap information) - this logic remains the same
for i, bar in enumerate(bars):
    width = bar.get_width()
    label = plot_data['Overlap'].iloc[i]
    ax.text(
        width + 0.1,
        bar.get_y() + bar.get_height() / 2,
        label,
        ha='left',
        va='center',
        fontsize=10,
        color=text_color
    )

# --- 4. Style the Plot for Publication Quality ---
# This section is identical, ensuring the same clean, academic style.

ax.set_title(
    'Top Enriched GO Biological Processes',
    fontsize=16, fontweight='bold', pad=20, color=text_color
)
ax.set_xlabel(
    '-log$_{10}$(Adjusted P-value)',
    fontsize=12, fontweight='medium', color=text_color
)
ax.set_ylabel('')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_color(grid_color)

ax.tick_params(axis='y', length=0, pad=10)
plt.xticks(fontsize=10, color=text_color)
plt.yticks(fontsize=11, color=text_color)
ax.grid(axis='x', linestyle='--', color=grid_color, alpha=0.7)
ax.set_axisbelow(True)
fig.tight_layout()

# --- 5. Save the Figure ---
file_path = "/root/Cycle/HarmoCycle_v1/Figure234_h5ad/go_enrichment_plot.pdf"
plt.savefig(file_path, bbox_inches='tight', dpi=300)

print(f"Plot saved to {file_path}")
plt.show()
