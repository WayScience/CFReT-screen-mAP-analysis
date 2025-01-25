#!/usr/bin/env python

# # Applying metrics

# In[1]:


import pathlib

# imports
import sys

import pandas as pd
from pycytominer.cyto_utils import load_profiles

sys.path.append("../../utils")
from utils import data_utils

# In[2]:


data_dir = pathlib.Path("../data")
results_dir = pathlib.Path("./results").resolve(strict=True)

# Setting path
metadata_cluster_path = (results_dir / "cluster/metadata_w_clusters.csv").resolve(strict=True)

# setting single-cell profile paths
profile_paths = list(data_dir.glob("*sc_feature_selected.parquet"))
if len(profile_paths) == 0:
    raise FileNotFoundError("Profiles were not found at the given directory")


# In[3]:


# loading in the data
shared_features = data_utils.find_shared_features(profile_paths)

# loading all single-cell profiles and updating it with the shared features
loaded_profiles_df = []
for single_cell_path in profile_paths:
    # loading in single cell feature selected data
    single_cell_df = load_profiles(single_cell_path)

    # append the updated profiles to the loaded_profiles_df
    loaded_profiles_df.append(single_cell_df[shared_features])

# Concatenate all the profiles
all_profiles_df = pd.concat(loaded_profiles_df, axis=0)

print(all_profiles_df.shape)
all_profiles_df.head()


# In[4]:


# parameters
metadata_treatments = "Metadata_treatment"
profile = None
target_name = None


# split the metadata and morphology feature
meta_cols, feat_cols = data_utils.split_meta_and_features(all_profiles_df)

# check if the selected metadata column contains the metadata_treatment that represents the control
if metadata_treatments not in meta_cols:
    raise ValueError(
        f"{metadata_treatments} is a metadata column that does not exist"
    )

target_df = all_profiles_df.loc[all_profiles_df[metadata_treatments] == target_name]
treated_df = all_profiles_df.loc[all_profiles_df[metadata_treatments] != target_name]
