"""
In this module, we have functions that has to do with providing utilts for our
"""

import pathlib
from collections import defaultdict

import pandas as pd
from pycytominer import annotate
from pycytominer.cyto_utils.features import infer_cp_features


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

def annotate_batch(
    profiles: dict, metadata_dir: str | pathlib.Path, platemap_metadata_query: str, profile_metadata_query: str
) -> dict:
    """ Annotate batch profiles with platemap metadata.

    This function takes a nested dictionary of profiles, annotates them with metadata
    from corresponding platemap files, and returns the annotated profiles.

    Parameters
    ----------
    profiles : dict
        A nested dictionary where:
        - Keys (str): Batch IDs.
        - Values (dict): Dictionaries mapping platemap names to lists of profile DataFrames.
    metadata_dir : str | pathlib.Path
        Path that points to the metadata folder
    platemap_metadata_query : str
        Column name in the platemap file to use as a join key.
    profile_metadata_query : str
        Column name in the profile DataFrame to use as a join key.

    Returns
    -------
    dict
        A nested dictionary where:
        - Keys (str): Batch IDs.
        - Values (dict): Dictionaries mapping platemap names to lists of annotated DataFrames.
    """
    # Initialize a dictionary to store annotated profiles for all batches
    annotated_batch_profiles = defaultdict(lambda: None)

    for batch_id, platemap in profiles.items():

        # Initialize a dictionary for annotated profiles for a single batch
        annotated_batch = {}

        for platemap_name, profile_list in platemap.items():
            # Construct the file path for the platemap metadata and load
            plate_map_path = (metadata_dir / f"{platemap_name}.csv").resolve(
                strict=True
            )
            platemap_df = pd.read_csv(plate_map_path)

            # Annotate each profile in the profile list
            annotated_profiles = []
            for profile in profile_list:
                annotated_df = annotate(
                    profiles=profile,
                    platemap=platemap_df,
                    join_on=[platemap_metadata_query, profile_metadata_query],
                )
                annotated_profiles.append(annotated_df)

            # Store the annotated profiles under the platemap name
            annotated_batch[platemap_name] = annotated_profiles

        # Store the batch-level annotations under the batch ID
        annotated_batch_profiles[batch_id] = annotated_batch

    return annotated_batch_profiles
