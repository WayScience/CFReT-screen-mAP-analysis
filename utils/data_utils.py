"""
In this module, we have functions that has to do with providing utilities for
the loaded profiles.
"""

import pathlib

import pandas as pd
import pyarrow.parquet as pq
from pycytominer.cyto_utils import infer_cp_features


def split_meta_and_features(
    profile: pd.DataFrame,
    compartments: list[str] = ["Nuclei", "Cells", "Cytoplasm"],
    metadata_tag: bool | None = False,
) -> tuple[list[str], list[str]]:
    """Splits metadata and feature column names

    Parameters
    ----------
    profile : pd.DataFrame
        Dataframe containing image-based profile
    compartments : list, optional
        compartments used to generated image-based profiles, by default
        ["Nuclei", "Cells", "Cytoplasm"]
    metadata_tag : Optional[bool], optional
        indicating if the profiles have metadata columns tagged with 'Metadata_'
        , by default False

    Returns
    -------
    tuple[List[str], List[str]]
        Tuple containing metadata and feature column names
    """

    # identify features names
    features_cols = infer_cp_features(profile, compartments=compartments)

    # iteratively search metadata features and retain order if the Metadata tag is not added
    if metadata_tag is False:
        meta_cols = [
            colname
            for colname in profile.columns.tolist()
            if colname not in features_cols
        ]
    else:
        meta_cols = infer_cp_features(profile, metadata=metadata_tag)

    return (meta_cols, features_cols)

def find_shared_features(profile_paths: list[str | pathlib.Path]) -> list[str]:
    """Find the shared features (columns) between the profiles in the provided list of
    file paths, while retaining the order of features as they appear in the first
    profile.

    This function leverages the schema information from the Parquet files to extract the
    the columns names without loading in the entire dataset. The first profile is used
    as the reference for the order of features. Then, the function iterates through the
    remaining profiles to find the common features. If no common features are found, an
    empty list is returned.

    Parameters
    ----------
    profile_paths : list[str | pathlib.Path]
        A list of file paths pointing to the Parquet profile files.

    Returns
    -------
    list[str]
        A list of features (column names) that are common across all profiles, retaining the order
        from the first profile.

    Raises
    ------
    ValueError
        If any of the profile paths do not point to Parquet files.
    """
    # type checker to check if the file provided are parquet files
    for path in profile_paths:
        if not path.suffix == ".parquet" or path.suffix == ".pq":
            raise ValueError("All profile paths must point to Parquet files.")

    # initialize the shared features to None
    shared_features = None

    # iterate through the profile paths
    for profile_path in profile_paths:
        # Load the metadata of the Parquet file
        parquet_metadata = pq.ParquetFile(profile_path)

        # Extract column names from the schema
        column_names = parquet_metadata.schema.names

        if shared_features is None:
            # Initialize shared features on the first iteration
            shared_features = column_names
        else:
            # Retain only the features that are shared, keeping their order
            shared_features = [name for name in shared_features if name in column_names]

    return shared_features if shared_features else []
