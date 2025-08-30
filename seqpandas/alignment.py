"""
Alignment visualization module for SAM/BAM files.
Shows sequences with alignment gaps relative to reference.
"""

import string
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from pysam import AlignmentFile

from .tools.pathing import pathing


def read_aligned_sam(
    alignment_file: str,
    reference_file: Optional[str] = None,
    region: Optional[str] = None,
    max_reads: Optional[int] = None,
    show_insertions: bool = True,
    gap_char: str = "-",
    insertion_char: str = "*",
) -> pd.DataFrame:
    """
    Read SAM/BAM file and return sequences aligned to reference coordinates.
    Each sequence will have the same length with gaps shown explicitly.

    Parameters
    ----------
    alignment_file : str
        Path to SAM/BAM file
    reference_file : str, optional
        Path to reference FASTA (for getting reference sequence)
    region : str, optional
        Region to fetch (e.g., "chr1:1000-2000")
    max_reads : int, optional
        Maximum number of reads to return
    show_insertions : bool
        If True, mark insertion positions; if False, skip insertions
    gap_char : str
        Character to use for deletions/gaps (default: '-')
    insertion_char : str
        Character to mark insertion positions (default: '*')

    Returns
    -------
    pd.DataFrame
        DataFrame with aligned sequences and metadata
    """
    alignment_path = pathing(alignment_file)

    # Determine file type and open
    if str(alignment_path).endswith(".bam"):
        samfile = AlignmentFile(str(alignment_path), "rb")
    else:
        samfile = AlignmentFile(str(alignment_path), "r")

    aligned_seqs = []

    # Get all alignments or just a region
    if region:
        alignments = samfile.fetch(region=region)
    else:
        alignments = samfile.fetch()

    # Process alignments
    for i, alignment in enumerate(alignments):
        if max_reads and i >= max_reads:
            break

        if alignment.is_unmapped:
            continue

        # Get aligned pairs (query_pos, ref_pos)
        aligned_pairs = alignment.get_aligned_pairs()

        # Build aligned sequence
        aligned_seq = build_aligned_sequence(
            alignment,
            aligned_pairs,
            show_insertions=show_insertions,
            gap_char=gap_char,
            insertion_char=insertion_char,
        )

        # Collect alignment info
        aligned_seqs.append(
            {
                "read_name": alignment.query_name,
                "ref_name": alignment.reference_name,
                "ref_start": alignment.reference_start,
                "ref_end": alignment.reference_end,
                "strand": "-" if alignment.is_reverse else "+",
                "mapq": alignment.mapping_quality,
                "aligned_sequence": aligned_seq,
                "original_sequence": alignment.query_sequence,
                "cigar": alignment.cigarstring,
            }
        )

    samfile.close()

    if not aligned_seqs:
        return pd.DataFrame()

    df = pd.DataFrame(aligned_seqs)

    # Pad sequences to same length (reference length in region)
    if len(df) > 0:
        df = pad_aligned_sequences(df, gap_char)

    return df


def build_aligned_sequence(
    alignment,
    aligned_pairs: List[tuple],
    show_insertions: bool = True,
    gap_char: str = "-",
    insertion_char: str = "*",
) -> str:
    """
    Build aligned sequence from alignment pairs.

    Parameters
    ----------
    alignment : pysam.AlignedSegment
        The alignment object
    aligned_pairs : List[tuple]
        List of (query_pos, ref_pos) pairs
    show_insertions : bool
        Whether to show insertions
    gap_char : str
        Character for gaps/deletions
    insertion_char : str
        Character for marking insertions

    Returns
    -------
    str
        Aligned sequence with gaps
    """
    if not alignment.query_sequence:
        return ""

    # Create position map
    ref_to_base = {}
    insertions = {}

    for query_pos, ref_pos in aligned_pairs:
        if ref_pos is not None and query_pos is not None:
            # Match or mismatch - check bounds
            if 0 <= query_pos < len(alignment.query_sequence):
                ref_to_base[ref_pos] = alignment.query_sequence[query_pos]
        elif ref_pos is None and query_pos is not None:
            # Insertion (query has base, ref doesn't)
            if show_insertions and aligned_pairs and 0 <= query_pos < len(alignment.query_sequence):
                # Find the last reference position
                last_ref = None
                for q, r in aligned_pairs[: aligned_pairs.index((query_pos, ref_pos))]:
                    if r is not None:
                        last_ref = r
                if last_ref is not None:
                    if last_ref not in insertions:
                        insertions[last_ref] = []
                    insertions[last_ref].append(alignment.query_sequence[query_pos])
        # If query_pos is None and ref_pos is not None, it's a deletion (gap in query)

    # Build the aligned sequence
    if not ref_to_base:
        return ""

    min_ref = min(ref_to_base.keys())
    max_ref = max(ref_to_base.keys())

    aligned = []
    for ref_pos in range(min_ref, max_ref + 1):
        if ref_pos in ref_to_base:
            aligned.append(ref_to_base[ref_pos])
            # Add insertions after this position
            if show_insertions and ref_pos in insertions:
                # Mark insertions with lowercase or special character
                for ins_base in insertions[ref_pos]:
                    aligned.append(ins_base.lower())
        else:
            # Deletion in query relative to reference
            aligned.append(gap_char)

    return "".join(aligned)


def pad_aligned_sequences(df: pd.DataFrame, gap_char: str = "-") -> pd.DataFrame:
    """
    Pad all aligned sequences to the same length based on reference coordinates.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with aligned sequences
    gap_char : str
        Character to use for padding

    Returns
    -------
    pd.DataFrame
        DataFrame with padded sequences
    """
    if len(df) == 0:
        return df

    # Find global start and end positions
    global_start = df["ref_start"].min()
    global_end = df["ref_end"].max()

    # Pad each sequence
    padded_seqs = []
    for _, row in df.iterrows():
        left_pad = gap_char * (row["ref_start"] - global_start)
        right_pad = gap_char * (global_end - row["ref_end"])
        padded_seq = left_pad + row["aligned_sequence"] + right_pad
        padded_seqs.append(padded_seq)

    df["aligned_sequence"] = padded_seqs
    df["alignment_start"] = global_start
    df["alignment_end"] = global_end

    return df


def show_alignment_block(
    df: pd.DataFrame, max_display: int = 10, width: int = 100, show_positions: bool = True
) -> str:
    """
    Create a text visualization of aligned sequences.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame from read_aligned_sam
    max_display : int
        Maximum number of sequences to show
    width : int
        Width of each line
    show_positions : bool
        Whether to show position numbers

    Returns
    -------
    str
        Formatted alignment block
    """
    if len(df) == 0:
        return "No alignments found"

    lines = []

    # Get sequences to display
    seqs_to_show = min(max_display, len(df))

    # Get alignment info
    if "alignment_start" in df.columns:
        start_pos = df["alignment_start"].iloc[0]
    else:
        start_pos = df["ref_start"].min()

    # Show position ruler
    if show_positions:
        lines.append(f"Position {start_pos}:")
        lines.append("=" * width)

    # Show each sequence
    for i in range(seqs_to_show):
        row = df.iloc[i]
        seq = row["aligned_sequence"]
        name = row["read_name"][:20] if len(row["read_name"]) > 20 else row["read_name"]

        # Show sequence in chunks
        for j in range(0, len(seq), width):
            chunk = seq[j : j + width]
            if j == 0:
                lines.append(f"{name:<22} {chunk}")
            else:
                lines.append(f"{'':<22} {chunk}")

    if len(df) > max_display:
        lines.append(f"\n... and {len(df) - max_display} more alignments")

    return "\n".join(lines)


# Integration with main seqpandas
def read_sam_aligned(
    alignment_file: str, reference_file: Optional[str] = None, **kwargs
) -> pd.DataFrame:
    """
    Convenience function for reading SAM/BAM with aligned view.

    Parameters
    ----------
    alignment_file : str
        Path to SAM/BAM file
    reference_file : str, optional
        Path to reference FASTA
    **kwargs
        Additional arguments for read_aligned_sam

    Returns
    -------
    pd.DataFrame
        DataFrame with aligned sequences
    """
    return read_aligned_sam(alignment_file, reference_file, **kwargs)


def alignment_to_frequency_matrix(
    alignment: List[str], normalize: bool = True, include_gaps: bool = True
) -> pd.DataFrame:
    """
    Convert aligned sequences to a frequency matrix.

    Parameters
    ----------
    alignment : List[str]
        List of aligned sequences (strings of equal length)
    normalize : bool
        If True, normalize counts to frequencies (0-1)
    include_gaps : bool
        If True, include gap characters in the frequency matrix

    Returns
    -------
    pd.DataFrame
        Frequency matrix with positions as rows and characters as columns
    """
    if not alignment:
        return pd.DataFrame()

    # Convert to uppercase for consistency
    alignment = [seq.upper() for seq in alignment]

    # Create DataFrame from alignment
    df = pd.DataFrame([list(seq) for seq in alignment])

    # Calculate frequency counts at each position
    counts_df = df.apply(lambda x: x.value_counts(), axis=0).fillna(0)

    # Get all unique characters
    characters = sorted(counts_df.index)

    # Ensure gaps are included if requested
    if include_gaps:
        if "-" not in characters:
            characters.append("-")
        if "." not in characters:
            characters.append(".")

    # Reindex to include all characters
    counts_df = counts_df.reindex(characters, fill_value=0)

    # Normalize if requested
    if normalize:
        freq_df = counts_df / counts_df.sum()
    else:
        freq_df = counts_df

    # Transpose for standard format (positions as rows)
    freq_df = freq_df.transpose()

    return freq_df


def generate_position_labels(
    reference_seq: str, positions: Optional[List[int]] = None
) -> Tuple[List[str], List[int]]:
    """
    Generate position labels for alignment with gap suffixes.

    Parameters
    ----------
    reference_seq : str
        Reference sequence (first sequence in alignment)
    positions : List[int], optional
        Specific positions to include (1-based)

    Returns
    -------
    Tuple[List[str], List[int]]
        Position labels and their indices
    """
    reference_seq = reference_seq.upper()
    ref_positions = []
    ref_counter = 1
    gap_counter = 0

    for i, char in enumerate(reference_seq):
        prefix = char
        next_char_gap = i + 1 < len(reference_seq) and reference_seq[i + 1] == "-"

        if char != "-":
            if next_char_gap:
                ref_positions.append(f"{prefix}{ref_counter}_{string.ascii_uppercase[0]}")
                gap_counter = 1
            else:
                ref_positions.append(f"{prefix}{ref_counter}")
            ref_counter += 1
        else:
            try:
                suffix = string.ascii_uppercase[gap_counter]
            except IndexError:
                suffix = "*"
                break
            ref_positions.append(f"{prefix}{ref_counter - 1}_{suffix}")
            gap_counter += 1

    # Filter positions if specified
    if positions:
        keep_indexes = []
        for i, element in enumerate(ref_positions):
            # Extract numeric part after optional letter prefix
            import re

            match = re.match(r"^[A-Za-z]?(\d+)", element)
            if match:
                numeric_value = int(match.group(1))
                if numeric_value in positions:
                    keep_indexes.append(i)

        filtered_positions = [ref_positions[i] for i in keep_indexes]
        return filtered_positions, keep_indexes

    return ref_positions, list(range(len(ref_positions)))


def get_alignment_frequencies(
    alignment_file: str,
    reference_file: Optional[str] = None,
    region: Optional[str] = None,
    max_reads: Optional[int] = None,
    positions: Optional[List[int]] = None,
    normalize: bool = True,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Read alignment file and generate frequency matrix with position labels.

    Parameters
    ----------
    alignment_file : str
        Path to SAM/BAM file
    reference_file : str, optional
        Path to reference FASTA
    region : str, optional
        Region to fetch
    max_reads : int, optional
        Maximum number of reads
    positions : List[int], optional
        Specific positions to include
    normalize : bool
        Normalize to frequencies

    Returns
    -------
    Tuple[pd.DataFrame, List[str]]
        Frequency DataFrame and position labels
    """
    # Read aligned sequences
    aligned_df = read_aligned_sam(
        alignment_file, reference_file=reference_file, region=region, max_reads=max_reads
    )

    if aligned_df.empty:
        return pd.DataFrame(), []

    # Get aligned sequences
    sequences = aligned_df["aligned_sequence"].tolist()

    # Generate frequency matrix
    freq_df = alignment_to_frequency_matrix(sequences, normalize=normalize)

    # Generate position labels using first sequence as reference
    ref_positions, indices = generate_position_labels(sequences[0], positions)

    # Filter frequency matrix if positions specified
    if positions and indices:
        freq_df = freq_df.iloc[indices, :]

    # Set column names to position labels
    freq_df.index = ref_positions

    return freq_df, ref_positions


def sequences_to_logo_format(
    sequences: List[str], reference_index: int = 0, positions: Optional[List[int]] = None
) -> Tuple[pd.DataFrame, List[str], List[int]]:
    """
    Convert sequences to format ready for logo plotting.

    Parameters
    ----------
    sequences : List[str]
        List of aligned sequences
    reference_index : int
        Index of reference sequence (default: 0)
    positions : List[int], optional
        Positions to include

    Returns
    -------
    Tuple[pd.DataFrame, List[str], List[int]]
        Frequency DataFrame transposed for logomaker, position labels, and variant positions
    """
    if not sequences:
        return pd.DataFrame(), [], []

    # Get frequency matrix
    freq_df = alignment_to_frequency_matrix(sequences)

    # Generate position labels
    ref_positions, indices = generate_position_labels(sequences[reference_index], positions)

    # Filter if positions specified
    if positions and indices:
        freq_df = freq_df.iloc[indices, :]

    # Transpose for logomaker format
    freq_df = freq_df.T
    freq_df.columns = ref_positions
    freq_df = freq_df[ref_positions]  # Ensure column order
    freq_df = freq_df.T.reset_index(drop=True)

    # Identify variant positions (where top frequency doesn't match reference)
    max_freq_cols = freq_df.idxmax(axis=1)
    variant_positions = []
    for i, (top_char, ref_pos) in enumerate(zip(max_freq_cols, ref_positions)):
        ref_char = ref_pos[0]  # First character is the reference base
        if top_char != ref_char:
            variant_positions.append(i)

    return freq_df, ref_positions, variant_positions
