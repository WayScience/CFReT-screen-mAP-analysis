"""
This module contains all the functions responsible for performing the analysis in this
project.
"""

import pathlib

import pandas as pd
from copairs import map

from .data_utils import load_config, split_meta_and_features


def calculate_dmso_batch_profiles():
    """"""


def calculate_trt_map_batch_profiles(
    batched_profiles: dict,
    configs: str | pathlib.Path | dict,
    outdir_path: str | pathlib.Path,
    shuffled: bool | None = False,
) -> None | tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate treatment-specific mean average precision (mAP) and average precision (AP) scores
    for batched phenotypic profiles.

    This function processes batches of phenotypic profiles, calculates average precision (AP)
    and mean average precision (mAP) scores, and saves the results to the specified output directory.
    It supports multiple control conditions to evaluate the phenotypic differences between treated
    and control samples.

    Parameters
    ----------
    batched_profiles : dict
        A dictionary where keys are batch identifiers and values are `pd.DataFrame` objects
        containing phenotypic profiles with metadata and features.

    configs : str, pathlib.Path, or dict
        Configuration settings for the analysis. If a path is provided, the configuration is
        loaded using the `load_config` function.

    outdir_path : str, pathlib.Path
        The directory where the AP and mAP scores will be saved. If None, no files are saved.

    shuffled : bool, optional
        Indicates whether the data should be shuffled. If set to True, a 'shuffled' label
        will be appended to the output file names. Default is False and provides 'original'
        label.

    Returns
    -------
    None
        AP and mAP scores are saved in provided directory path.

    Raises
    ------
    TypeError
        If the input types of `batched_profiles`, `configs`, or `outdir_path` are invalid.
    """

    # type checking
    if not isinstance(batched_profiles, dict):
        raise TypeError("'batched_profiles' must be a dictionary")
    # load configs if a path is provided and then load it
    if isinstance(configs, (str, pathlib.Path)):
        configs = load_config(configs)
    if not isinstance(configs, dict):
        raise TypeError("'configs' must be a dictionary")
    if not isinstance(outdir_path, (str, pathlib.Path)):
        raise TypeError("'outdir_path' must be either pathlib.Path or str")
    if not isinstance(shuffled, bool):
        raise TypeError("'shuffled' must be a boolean")

    # setting in project configs
    general_configs = configs["general_configs"]
    copairs_ap_configs = configs["copairs_ap_configs"]
    copairs_map_configs = configs["copairs_map_configs"]

    # setting shuffled labels
    shuffled_label = "original"
    if shuffled:
        shuffled_label = "shuffled"

    # Define control conditions for the analysis
    # Each tuple specifies the control type, treatment, and associated cell state
    control_list = [("negative", "DMSO", "failing"), ("positive", "DMSO", "healthy")]

    # Iterate over each batch of loaded plate profiles
    for batch_id, profile in batched_profiles.items():
        # Analyze the profile for each control condition
        for control_type, control_treatment, cell_state in control_list:
            # Create a copy of the profile to preserve the original data
            profile = profile.copy()

            # Assign a default reference index based on the row index
            profile["Metadata_reference_index"] = profile.index

            # Mark all non-control replicates (e.g., treatments not matching the current control)
            profile.loc[
                (profile["Metadata_treatment"] != control_treatment)
                & (profile["Metadata_cell_type"] != cell_state),
                "Metadata_reference_index",
            ] = -1

            # Move the "Metadata_reference_index" column to the beginning for clarity
            profile.insert(
                0, "Metadata_reference_index", profile.pop("Metadata_reference_index")
            )

            # Separate metadata columns from feature columns for downstream calculations
            meta_columns, feature_columns = split_meta_and_features(profile)

            # Calculate average precision (AP) for the profile
            # Positive pairs are based on treatments with the same metadata
            # Negative pairs compare all DMSO-treated wells to all treatments
            replicate_aps = map.average_precision(
                meta=profile[meta_columns],
                feats=profile[feature_columns].values,
                pos_sameby=copairs_ap_configs["pos_sameby"],
                pos_diffby=copairs_ap_configs["pos_diffby"],
                neg_sameby=[],
                neg_diffby=copairs_ap_configs["neg_diffby"],
            )

            # Exclude wells treated with the control treatment (DMSO)
            replicate_aps = replicate_aps.loc[
                replicate_aps["Metadata_treatment"] != control_treatment
            ]

            # Save the calculated AP scores to a file for further analysis
            save_ap_path = (
                outdir_path
                / f"{batch_id}_{shuffled_label}_{control_type}_control_{cell_state}_{control_treatment}_AP_scores.csv"
            )
            replicate_aps.to_csv(
                save_ap_path,
                index=False,
            )

            # Calculate mean average precision (mAP) from the AP scores
            replicate_maps = map.mean_average_precision(
                replicate_aps,
                sameby=copairs_map_configs["same_by"],  # Grouping criteria for mAP
                null_size=copairs_map_configs["null_size"],  # Null distribution size
                threshold=copairs_map_configs["threshold"],  # Significance threshold
                seed=general_configs["seed"],  # Seed for reproducibility
            )

            # Save the mAP scores to a file for reporting
            save_map_path = (
                outdir_path
                / f"{batch_id}_{shuffled_label}_{control_type}_control_{cell_state}_{control_treatment}_mAP_scores.csv"
            )
            replicate_maps.to_csv(
                save_map_path,
                index=False,
            )
