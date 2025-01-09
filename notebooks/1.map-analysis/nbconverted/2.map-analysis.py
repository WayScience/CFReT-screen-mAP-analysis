#!/usr/bin/env python

# ## 2. mAP analysis
#
# In this notebook we are analyzing the mAP scores generate from the previous step shown [here](./1.run-map.ipynb). In this note book, we are investigating that mAP scores generated to see if we can identify specific compounds that show evidence of reversing cardiac fibroblast.

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Setting up input and output paths

# In[2]:


results_dir = pathlib.Path("./results").resolve(strict=True)
map_results_dir = (results_dir / "map_scores").resolve(strict=True)

# Setting AP and mAP scores paths
nc_dmso_ap_path = (map_results_dir / "negative_control_failing_DMSO_AP_scores.csv").resolve(strict=True)
pc_dmso_ap_path = (map_results_dir / "positive_control_healthy_DMSO_AP_scores.csv").resolve(strict=True)
nc_dmso_map_path = (map_results_dir / "negative_control_failing_DMSO_mAP_scores.csv").resolve(strict=True)
pc_dmso_map_path = (map_results_dir / "positive_control_healthy_DMSO_mAP_scores.csv").resolve(strict=True)

# make directory for figures
map_analysis_results_dir = (results_dir / "map_analysis").resolve()
map_analysis_results_dir.mkdir(exist_ok=True)

fig_dir_path = (map_analysis_results_dir / "figures").resolve()
fig_dir_path.mkdir(exist_ok=True)

# pathway informaiton
pathway_info = pathlib.Path("../data/metadata/original_platemaps/pathways_platemap.csv").resolve(strict=True)


# Here, we are loading the mAP scores generated using both positive and negative controls as references. These two datasets are then merged into a single DataFrame to combine the scores for each treatment. Next, we calculate the delta mAP, which is the difference between the mAP score generated using the negative control (failing CF cells + DMSO) as the reference and the mAP score generated using the positive control (healthy CF cells + DMSO) as the reference.

# In[3]:


interested_cols = ["Metadata_treatment", "mean_average_precision", "corrected_p_value"]
pathway_df = pd.read_csv(pathway_info)[["UCD ID", "Pathway"]]
nc_map_df = pd.read_csv(nc_dmso_map_path)[interested_cols]
pc_map_df = pd.read_csv(pc_dmso_map_path)[interested_cols]

# rename columns
nc_map_df = nc_map_df.rename(columns={"mean_average_precision":"negative_mean_average_precision", "corrected_p_value":"negative_corrected_p_value"})
pc_map_df = pc_map_df.rename(columns={"mean_average_precision":"positive_mean_average_precision", "corrected_p_value":"positive_corrected_p_value"})

# merge both together
all_map_df = nc_map_df.merge(pc_map_df, on="Metadata_treatment", how="inner")

# calculating delta mAP
all_map_df["delta_mAP"] = (
    all_map_df["negative_mean_average_precision"]
    - all_map_df["positive_mean_average_precision"]
)

# now leta add pathway information into this data
all_map_df = all_map_df.merge(pathway_df, left_on="Metadata_treatment", right_on="UCD ID", how="inner")
all_map_df = all_map_df.drop(columns="UCD ID")
all_map_df.head()


# ## Dual-Control mAP Analysis: Evaluating Therapeutic Potential Through Comparative Phenotypic Scoring

# In this section, we visualize both mAP scores using a scatter plot to examine the distribution of points and interpret the relationships between the scores generated from the negative and positive controls. This helps us explore how treatments compare across these two references and identify trends or patterns in the data.

# In[4]:


# Create a jointplot
dual_map_plot = sns.jointplot(
    data=all_map_df,
    x="negative_mean_average_precision",
    y="positive_mean_average_precision",
    hue="Pathway",
    palette="tab10",
    kind="scatter",
    marginal_kws=dict(fill=True, alpha=0.4),
    joint_kws=dict(s=110, alpha=0.8, edgecolor="black"),
)
# Make the figure larger
dual_map_plot.fig.set_size_inches(10, 8)

# Add a diagonal line
dual_map_plot.ax_joint.plot([0, 1], [0, 1], color="red", linestyle="--", linewidth=1.5)

# Set axis limits
dual_map_plot.ax_joint.set_xlim(0, 1.02)
dual_map_plot.ax_joint.set_ylim(0, 1.02)

# Add gridlines to the central plot
dual_map_plot.ax_joint.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.7)

# Add axis labels
dual_map_plot.ax_joint.set_xlabel("mAP Score (Reference: Failing CF Cells + DMSO)", fontsize=12)
dual_map_plot.ax_joint.set_ylabel("mAP Score (Reference: Healthy CF Cells +DMSO)", fontsize=12)

# Add the overall title
dual_map_plot.fig.suptitle("mAP Scores: Negative vs. Positive Controls", fontsize=16, y=1.0)

# Adjust layout to avoid overlap
dual_map_plot.fig.tight_layout()

# Adjust legend size
legend = dual_map_plot.ax_joint.legend_
legend.set_title("Pathway")
legend.set_bbox_to_anchor((1.51, 1))  # Position the legend outside the plot (optional)
for text in legend.get_texts():
    text.set_fontsize(10)  # Set font size for legend labels
legend.get_title().set_fontsize(12)  # Set font size for legend title

# Save the figure
plt.savefig(fig_dir_path / "map_scores_with_density.png", dpi=300, bbox_inches="tight")

# Show the plot
plt.show()


# ## Differential mAP scores Histogram analysis (Pathway Level)
#
# The primary goal of using differential mAP scores (delta_MAP) in this high-content screening is to evaluate the therapeutic potential of 550 compounds for reversing cardiac fibrosis. By comparing the morphological profiles of treated diseased cells against two controls failing cardiac cells + DMSO (negative control) and healthy cardiac cells + DMSO (positive control) the differential mAP scores help identify compounds that effectively shift the diseased phenotype toward the healthy state.

# $$
# \begin{equation}
# \Delta \text{mAP} = \text{mAP}_{\text{Negative Control}} - \text{mAP}_{\text{Positive Control}}
# \end{equation}
# $$
#
# $\text{mAP}_{\text{Negative Control}}  ( X )$:
#
# - The mean average precision score when comparing treated diseased cells to the failing control.
# - A **lower** $X$ indicates that the compound’s effect is similar to the diseased phenotype (minimal therapeutic relevance).
# - A **higher**$X$ indicates that the compound’s effect is different to the diseased phenotype (potential therapeutic relevance)
#
# $\text{mAP}_{\text{Positive Control}}  ( Y )$:
#
# - The mean average precision score when comparing treated diseased cells to the healthy control.
# - A **lower** $Y$ indicates that the compound shifts the disease state closer to to the healthy type (potential therapeutic relevance)
# - A **higher** $Y$ indicates that the compound shifts the diseased cells further to the healthy phenotype. (minimal therapeutic relevance)
#
# $ \Delta \text{mAP}  ( X - Y ):$
# - The difference between  $\text{mAP}{\text{Positive Control}}  and  \text{mAP}{\text{Negative Control}}$ .
# - This score provides a measure of how effectively a compound reverses the failing phenotype in a single value.
#
# #### What is the range of the scores?
#
# The range of these scores is [-1, 1]
#
# - **Positive** $\Delta \text{mAP}$ $(0 < \Delta \text{mAP} \leq 1)$:
#   Indicates the compound aligns more with the healthy phenotype than the diseased phenotype.
#
# - **Negative** $\Delta \text{mAP}$ $(-1 \leq \Delta \text{mAP} < 0)$:
#   Indicates the compound aligns more with the diseased phenotype than the healthy phenotype.
#
# - **Zero** $\Delta \text{mAP}$ $(\Delta \text{mAP} = 0)$:
#   Indicates no difference between the compound's effect on healthy and diseased phenotypes.
#

# In[5]:


# Plot stacked bar histogram using seaborn
plt.figure(figsize=(10, 6))
sns.histplot(
    data=all_map_df,
    x="delta_mAP",
    hue="Pathway",
    multiple="stack",
    bins=20,
    palette="tab10",
    alpha=0.8
)

# Add plot labels and title
plt.title("Delta mAP Distribution by Pathway", fontsize=16)
plt.xlabel("Delta mAP (Positive mAP - Negative mAP)", fontsize=12)
plt.ylabel("Frequency", fontsize=12)
plt.tight_layout()

# save plot
plt.savefig(fig_dir_path / "delta_mAP_histogram.png", dpi=300, bbox_inches="tight")

# Show the plot
plt.show()


# Leveraging Delta mAP can also provide an opportunity to rank them from highest to lowest and identify potential hits.

# In[6]:


plt.figure(dpi=200)

# only getting the delta map scores
delta_map_df = all_map_df[["Metadata_treatment", "delta_mAP", "Pathway"]]

# Add ranks to the DataFrame
delta_map_df["rank"] = delta_map_df["delta_mAP"].rank(ascending=True, method="max")
delta_map_df = delta_map_df.sort_values(by="rank", ascending=False)

sns.scatterplot(
    data=delta_map_df, x="rank", y="delta_mAP", hue="Pathway", palette="tab10"
)

axis_padding = 0.05
plt.title("delt_mAP Rankings", fontsize=14)
plt.xlim(0, delta_map_df["rank"].max() + 1)
plt.ylim(-axis_padding + delta_map_df["delta_mAP"].min(), delta_map_df["delta_mAP"].max() + axis_padding)
plt.rc("legend", fontsize=5)
plt.rc("legend", title_fontsize=10)
plt.legend( loc="upper left")

# save ranked delta maps
delta_map_df.to_csv(map_analysis_results_dir / "ranked_delta_mAPs.csv", index=False)

# save plot
plt.savefig(fig_dir_path / "delta_mAP_rankings.png", dpi=300, bbox_inches="tight")

# display plot
plt.show()
