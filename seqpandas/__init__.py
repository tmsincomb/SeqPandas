"""Top-level package for SeqPandas."""

from .alignment import (
    alignment_to_frequency_matrix,
    generate_position_labels,
    get_alignment_frequencies,
    read_aligned_sam,
    read_sam_aligned,
    sequences_to_logo_format,
    show_alignment_block,
)
from .bed import read_bed
from .seqpandas import BioDataFrame, read, read_seq
from .tools.pathing import pathing
from .vcf import read_vcf

# Optional pileup module - requires scipy
try:
    from .pileup import Pileup

    __all__ = [
        "read",
        "read_seq",
        "BioDataFrame",
        "read_vcf",
        "read_bed",
        "pathing",
        "read_aligned_sam",
        "read_sam_aligned",
        "show_alignment_block",
        "alignment_to_frequency_matrix",
        "generate_position_labels",
        "get_alignment_frequencies",
        "sequences_to_logo_format",
        "Pileup",
    ]
except ImportError:
    __all__ = [
        "read",
        "read_seq",
        "BioDataFrame",
        "read_vcf",
        "read_bed",
        "pathing",
        "read_aligned_sam",
        "read_sam_aligned",
        "show_alignment_block",
        "alignment_to_frequency_matrix",
        "generate_position_labels",
        "get_alignment_frequencies",
        "sequences_to_logo_format",
    ]

__author__ = """Troy Sincomb"""
__email__ = "troysincomb@gmail.com"
__version__ = "0.0.2"
