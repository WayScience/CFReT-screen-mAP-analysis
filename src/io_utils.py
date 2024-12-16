"""io_utils.py

This module will contain functions pertaining to loading and writing files.
"""

import pathlib
from collections import defaultdict

import pandas as pd
import yaml


def load_config(fpath: str | pathlib.Path) -> dict:
    """Loads in configuration file if specified path

    Parameters
    ----------
    fpath : str | pathlib.Path
        path pointing to config file

    Returns
    -------
    dict
        contents of the configuration file
    """

    # type checking
    if not isinstance(fpath, (str | pathlib.Path)):
        raise TypeError("'fpath' must be a string or pathlib.Path object")
    if isinstance(fpath, str):
        fpath = pathlib.Path(fpath)

    # check if the path exists and returns the full path
    fpath = fpath.resolve(strict=True)

    # next is to load the yaml file
    with open(fpath) as content:
        return yaml.safe_load(content)


def load_barcodes(barcode_path: str | pathlib.Path) -> dict:
    """
    Load barcode data from a given file and organize it into batches.

    This function reads barcode data, groups it by the 'platemap_file' column, and
    organizes the grouped data into batches. Each batch is assigned a unique batch ID
    (e.g., 'batch_1', 'batch_2', etc.).

    Parameters
    ----------
    barcode_path : str | pathlib.Path
        Path to the file containing barcode data. This path can be a string or a pathlib.Path
        object.

    Returns
    -------
    dict
        A dictionary where keys are batch IDs (e.g., "batch_1", "batch_2") and values are
        dictionaries.  Each inner dictionary maps a 'platemap_file' name to the
        corresponding 'plate_barcode' data.
    """

    # type checking
    if not isinstance(barcode_path, (str | pathlib.Path)):
        raise TypeError("'barcode_path' must a str or pathlib.Path")
    if isinstance(barcode_path, str):
        barcode_path = pathlib.Path(barcode_path)

    # load barcode as csv
    barcodes = pd.read_csv(barcode_path)

    # Initialize an empty dictionary to store barcode contents
    barcode_contents = defaultdict(lambda: None)

    # Group barcodes by 'platemap_file' and iterate through each group
    for idx, (platemap_name, df) in enumerate(
        barcodes.groupby("platemap_file"), start=1
    ):
        # Generate a batch ID
        batch_id = f"batch_{idx}"

        # Collect plate barcodes for the current platemap
        batch_content = {platemap_name: df["plate_barcode"].values.tolist()}

        # Store the batch content in the main dictionary
        barcode_contents[batch_id] = batch_content

    return barcode_contents
