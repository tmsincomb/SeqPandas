"""Top-level package for SeqPandas."""
from .seqpandas import read, read_seq, BioDataFrame
from .vcf import read_vcf
from .bed import read_bed
from .pathing import pathing


__author__ = """Troy Sincomb"""
__email__ = "troysincomb@gmail.com"
__version__ = "0.0.2"
__all__ = ["read", "read_seq", "BioDataFrame", "read_vcf", "read_bed", "pathing"] 
