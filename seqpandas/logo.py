import argparse
import json
import logging
import multiprocessing
from pathlib import Path
from typing import List, Literal

import logomaker
import matplotlib.pyplot as plt
import pandas as pd
import patchworklib as pw
from matplotlib import ticker as mtick

# Import alignment functions
from .alignment import sequences_to_logo_format

# Import the custom glyphs from JSON

# Use __file__ to create a path relative to the current file
CUSTOM_GLYPHS_PATH = Path(__file__).parent / "colors" / "custom_glyphs.json"
with open(CUSTOM_GLYPHS_PATH) as f:
    CUSTOM_GLYPHS = json.load(f)

DEFAULT_THREADS = multiprocessing.cpu_count()
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def swap_elements(my_list, index1, index2):
    """
    Swap two elements in a list at the specified indices.

    Args:
        my_list (list): The list containing elements to swap
        index1 (int): Index of the first element
        index2 (int): Index of the second element

    Returns:
        list: The list with swapped elements

    Raises:
        ValueError: If either index is out of range
    """
    # Create a copy to avoid modifying the original
    result = my_list.copy()

    # Check if indices are valid
    if 0 <= index1 < len(result) and 0 <= index2 < len(result):
        # Swap elements at the specified indices
        result[index1], result[index2] = result[index2], result[index1]
        return result
    else:
        raise ValueError(f"Invalid indices: {index1}, {index2} for list of length {len(result)}")


def plot_logo_with_indels(
    alignment: list[str],
    output_file: Path | str,
    title: str = "Sequence Logo with Indels",
    positions: list[int] | None = None,
    highlight_pos: list[int] = [],
    highlight_color: str = "gold",
    highlight_alpha: float = 0.5,
    color_scheme: Literal[
        "NajafabadiEtAl2017", "PGDM1400", "PGT121", "VRC01", "VRC07.523LS", "grays"
    ] = "NajafabadiEtAl2017",
):
    """
    Create a sequence logo plot accounting for insertions (gaps '-') and deletions ('.').
    Uses the first sequence as the reference for x-ticks and displays alphabetical suffixes for positions
    followed by gaps or redundant indexes.

    Parameters:
    alignment (list of str): List of aligned sequences (strings of equal length).
    title (str): Title of the plot.
    """
    # check if positions is list of string, if so convert to int
    if positions:
        try:
            positions = [int(pos) for pos in positions]
        except ValueError:
            logger.info("Positions should be a list of integers.")
            raise ValueError("Positions should be a list of integers.")

    # Use the new alignment functions to get frequency matrix and labels
    freq_df, ref_positions, variant_positions = sequences_to_logo_format(
        alignment, reference_index=0, positions=positions
    )

    # Use variant positions as highlight positions if not specified
    if not highlight_pos:
        highlight_pos = variant_positions

    # Handle custom color schemes
    special_char_colors = []
    if color_scheme in ["PGDM1400", "PGT121", "VRC01", "VRC07.523LS"]:
        special_colors = CUSTOM_GLYPHS[color_scheme]
        for i, ref_pos in enumerate(ref_positions):
            ref_char = ref_pos.split("-")[0]
            if ref_char in special_colors:
                color = special_colors[ref_char]
                special_char_colors.append((ref_char[0], i, color))
        # make backgroud gray
        dark_gray_scheme = {
            "A": "#444444",
            "C": "#3A3A3A",
            "D": "#2F2F2F",
            "E": "#4A4A4A",
            "F": "#383838",
            "G": "#474747",
            "H": "#3E3E3E",
            "I": "#414141",
            "K": "#2B2B2B",
            "L": "#323232",
            "M": "#353535",
            "N": "#2D2D2D",
            "P": "#404040",
            "Q": "#303030",
            "R": "#272727",
            "S": "#363636",
            "T": "#3B3B3B",
            "V": "#333333",
            "W": "#292929",
            "Y": "#393939",
        }
        color_scheme = dark_gray_scheme

    # Save frequency matrix as CSV
    df = freq_df.copy().T
    df.columns = ref_positions
    df = round(df, 3)
    df.to_csv(Path(output_file).with_suffix(".csv"), index=True)

    # Create a visualization
    ax = pw.Brick(figsize=(max(8, len(ref_positions) // 2), 0.5))

    # breakpoint()
    ww_logo = logomaker.Logo(freq_df, ax=ax, color_scheme=color_scheme, show_spines=False)

    for c, p, color in special_char_colors:
        ww_logo.style_single_glyph(c=c, p=p, color=color)
        # ww_logo.style_single_glyph(c="K", p=0, color="blue")
    for pos in highlight_pos:
        ww_logo.highlight_position(pos, color=highlight_color, alpha=highlight_alpha)

    ax.xaxis.set_major_locator(mtick.FixedLocator(range(0, len(ref_positions))))
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Reference Position")
    ax.set_ylabel("Frequency")
    ax.set_xticklabels(labels=ref_positions)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.savefig(output_file, bbox_inches="tight", dpi=300)
    return freq_df


def process_file(
    input_df_path: str, output_logo: str, positions: list[int], ref_pos: int, color_scheme: str
) -> None:
    """Process a single DataFrame file and generate a logo."""
    logger.info(f"Processing file: {input_df_path}")
    aln_df = pd.read_csv(input_df_path, index_col=0)
    # Handle multi-character entries
    aln_df = aln_df.map(lambda x: "-" if pd.notna(x) and len(str(x)) > 1 else x)
    pileup = aln_df.apply(lambda x: "".join(x), axis=1).values
    # replace first index with ref_pos
    pileup = swap_elements(pileup, 0, ref_pos)
    plot_logo_with_indels(
        alignment=pileup, output_file=output_logo, positions=positions, color_scheme=color_scheme
    )


def process_folder(
    input_folder: str,
    output_folder: str,
    positions: list[int],
    file_pattern: str = "*.tsv",
    ref_pos: int = 0,
    color_scheme: str = "NajafabadiEtAl2017",
) -> None:
    """Process all DataFrame files in a folder and generate logos."""
    input_path = Path(input_folder)
    output_path = Path(output_folder)

    # Find all matching files in the input folder
    files = list(input_path.glob(file_pattern))

    if not files:
        logger.warning(f"No matching files found in {input_folder} using pattern {file_pattern}")
        return

    logger.info(f"Found {len(files)} files to process")

    for file_path in files:
        # Generate output filename
        output_file = output_path / f"{file_path.stem}_logo.png"
        process_file(str(file_path), str(output_file), positions, ref_pos, color_scheme)


def main():
    parser = argparse.ArgumentParser(description="Plot sequence logo from alignment DataFrame")
    parser.add_argument(
        "input",
        type=str,
        help="Input DataFrame file or folder containing DataFrame files",
    )
    parser.add_argument("output", type=str, help="Output logo file or folder for output logos")
    parser.add_argument("--positions", type=int, nargs="+", help="Positions to plot")
    parser.add_argument(
        "--file-pattern",
        type=str,
        default="*.csv",
        help="File pattern to match when processing folders (default: *.tsv)",
    )
    parser.add_argument(
        "--is-folder", action="store_true", help="Explicitly indicate input is a folder"
    )
    parser.add_argument(
        "--ref-pos",
        type=int,
        default=0,
        help="Where the reference is located in the aln file",
    )
    parser.add_argument(
        "--color-scheme",
        type=str,
        default="NajafabadiEtAl2017",
        help="Color theme for the logo (default: NajafabadiEtAl2017)",
    )
    args = parser.parse_args()

    # Determine if input is a file or folder
    input_path = Path(args.input)
    is_folder = args.is_folder or input_path.is_dir()

    if is_folder:
        logger.info(f"Processing folder: {args.input}")
        process_folder(
            args.input,
            args.output,
            args.positions,
            args.file_pattern,
            args.ref_pos,
            args.color_scheme,
        )
    else:
        logger.info(f"Processing file: {args.input}")
        process_file(args.input, args.output, args.positions, args.ref_pos)


if __name__ == "__main__":
    main()


# Example usage:
#
# Process a single file:
#   python logo.py input_alignment.tsv output_logo.png --positions 1 2 3 4 5
#
# Process all TSV files in a folder:
#   python logo.py input_folder/ output_folder/ --positions 1 2 3 4 5 --is-folder
#
# Process specific file pattern in a folder:
#   python logo.py input_folder/ output_folder/ --positions 1 2 3 4 5 --file-pattern "*_alignment.tsv" --is-folder
