"""Main module."""

import bz2
import gzip
from abc import ABC
from copy import deepcopy  # copies all nested dictionaries
from io import StringIO
from pathlib import Path
from types import GeneratorType
from typing import Iterator, List, Tuple, Union

import pandas as pd
from Bio import SeqIO  # __init__ kicks in
from Bio.Seq import Seq
from pysam import AlignmentFile

from .tools.pathing import pathing

# from pandas.core.frame import * # Needed if I want to dive even deeper.


class SubclassedSeries(pd.Series):
    """Pandas Series API to Inherit"""

    @property
    def _constructor(self):
        return SubclassedSeries

    @property
    def _constructor_expanddim(self):
        return SubclassedDataFrame


class SubclassedDataFrame(pd.DataFrame, ABC):
    """Pandas DataFrame to Inherit"""

    @property
    def _constructor(self):
        return SubclassedDataFrame

    @property
    def _constructor_sliced(self):
        return SubclassedSeries


# todo give a dataframe a method extension to become a tfrecord https://www.tensorflow.org/tutorials/load_data/tfrecord
class BioDataFrame(SubclassedDataFrame, ABC):
    """Expanded Pandas DataFrame to handle BioPython SeqRecords generator or genomic file types"""

    @classmethod
    def from_seqrecords(
        cls,
        seqrecords: Union[GeneratorType, list],
        *args,
        **kwargs,
    ) -> pd.DataFrame:
        """Takes Biopython parsed output to convert to a proper DataFrame
            See from_records for more details on the rest of the parameters
        :param seqrecords: Generator or list from BioPython universal output. All formats are the same output.
        >>> from_seqrecords(Bio.SeqIO.parse('file.fasta', format='fasta'))
        """
        # Nomalize nested featues; will result in some redundant.
        # This won't be an issue since small seqrecord count usually means high
        # amount of features & vice versa .
        data = cls.__normalize_seqrecords(seqrecords)
        return cls.from_records(data, *args, **kwargs)

    def to_vcf(self):
        return self

    @staticmethod
    def __normalize_seqrecords(seqrecords: Union[GeneratorType, list]) -> List[dict]:
        """Pull nested dictionaries into a single dictionary.
        Priority is given to the keys higher in the hierarchy.
        :param seqrecords: Generator from BioPython universal output. All formats are the same output.
        :returns: List of dictionaries with keys that were obtained along the way.
        >>> __normalize_seqrecords(Bio.SeqIO.parse('file.fasta', format='fasta'))
        """

        def update_record(_record, data, reference_keys):
            """Quick update of dictionary without updating prexisting keys"""
            for key, value in data.items():
                if key not in reference_keys:
                    _record[key] = value
            return _record

        # Handle SAM/BAM/CRAM alignment dicts directly
        seqrecords_list = list(seqrecords)
        if seqrecords_list and isinstance(seqrecords_list[0], dict):
            return seqrecords_list

        seqrecords = SeqIO.to_dict(
            seqrecords_list
        ).values()  # Only care for the actual records themselves

        records = []
        for seqrecord in seqrecords:
            _records = []
            # SeqIO Parse is a nested class
            record = seqrecord.__dict__
            # If a more complicated format is used; features will be nested.
            features = record.pop("features") if record.get("features") else []
            # Annotation dictionary inside each seqrecord
            annotations = record.pop("annotations") if record.get("annotations") else {}
            # Add each annotation as seperate column
            record = update_record(record, annotations, record.keys())
            for feature in features:
                _record = deepcopy(record)  # Needs to refresh for each feature
                # Meta that make up the feature
                aspects = feature.__dict__
                # Qualifier dictionary inside each feature
                qualifiers = aspects.pop("qualifiers") if aspects.get("qualifiers") else {}
                # Add each feature aspect
                _record = update_record(_record, aspects, record.keys())
                # Add each qualifier
                _record = update_record(_record, qualifiers, record.keys())
                # Collect normalized feature
                _records += [_record]
            # If no normalized feature collected use original seq record
            if not _records:
                _records += [record]
            # Add current records list to past iterations.
            # We do this because there could be more than one feature per seqrecord.
            records += _records

        return records


def open_mode(suffix: str) -> Tuple[str, str]:
    """
    Returns the file compression mode and open mode for a given file suffix.

    Parameters
    ----------
    suffix : str
        File suffix, but specifically the last suffix of the file name.

    Returns
    -------
    Tuple[str, str]
        File compression mode and open mode.
    """
    suffix2open_mode = {
        ".gz": (gzip.open, "rt"),
        ".bz2": (bz2.open, "rt"),
        ".ab1": (open, "rb"),
    }
    return suffix2open_mode.get(suffix, (open, "r"))


def pull_format(path: Path, override: str = None) -> str:
    """
    Pulls the Bioinformatics file format from the path suffix. This is for BioPython
    so it will be the suffix without the "." I.E. "fasta" or "fastq" and not ".fasta" or ".fastq"

    Parameters
    ----------
    path : Path
        Path to the file.
    override : str, optional
        Format is known so there is no need to pull last suffix from path, by default None

    Returns
    -------
    str
        Format of the file. I.E. fasta, fastq
    """
    short2format = {
        "fa": "fasta",
        "fq": "fastq",
    }  # short form suffixes commonly used in bioinformatics
    # Get file format with priority: input, file extension, or default fasta
    format = override or (path.suffixes[0] if path.suffix else None) or "fasta"
    # BioPython needs a perfect lowercase format name
    format = format.lower().strip().replace(".", "")
    # replace shorthand file types
    format = short2format.get(format, format)
    return format


def read(
    handle: Union[str, StringIO], format: str = None, simple: bool = False
) -> Iterator[Union[SeqIO.SeqRecord, tuple]]:
    """Reads a genomic file and converts it to SeqRecords.

    Parameters
    ----------
    handle : str or StringIO
        Path to a genomic file or a StringIO object.
    format : str, optional
        File format. Default is fasta
    simple : bool, optional
        If True, returns a generator of tuples instead of SeqRecords (no annotation objects). Giant speed boost!

    Returns
    -------
    pandas.DataFrame

    Examples
    --------
    >>> read('file.fasta.gz')
    >>> read('file.fasta.bz2')
    >>> read('file.fa.gz')
    >>> read('file.fq.gz')
    >>> read('nosuffix_file', 'fasta')
    """

    # Clean file path - make sure it exists
    if not isinstance(handle, StringIO):
        path = pathing(handle)
        format = pull_format(path, override=format)
    else:
        path = None
        if not format:
            raise ValueError("Must specify format for StringIO")
        format = format.lower().strip()

    # Special handling for SAM/BAM/CRAM files - they yield alignment dicts
    if format in ["sam", "bam", "cram"]:
        if format == "sam":
            samfile = AlignmentFile(str(path) if path else handle, "r", threads=1)
        elif format == "bam":
            samfile = AlignmentFile(str(path) if path else handle, "rb", threads=1)
        elif format == "cram":
            samfile = AlignmentFile(str(path) if path else handle, "rc", threads=1)

        # Yield alignments as dict-like objects
        for alignment in samfile:
            yield alignment.to_dict()
        samfile.close()
        return

    def read_handle(handle, format):
        if simple and format == "fastq":
            seqs: Iterator[tuple] = SeqIO.QualityIO.FastqGeneralIterator(handle)
            for annotation, sequence, quality_scores in seqs:
                yield annotation, Seq(sequence), quality_scores
        elif simple and format == "fasta":
            seqs: Iterator[tuple] = SeqIO.FastaIO.SimpleFastaParser(handle)
            for annotation, sequence, quality_scores in seqs:
                yield annotation, Seq(sequence), quality_scores
        else:
            seqs: SeqIO.SeqRecord = SeqIO.parse(handle, format=format)
            for seq in seqs:
                yield seq

    # Special case for StringIO if the file is already opened
    if isinstance(handle, StringIO):
        yield from read_handle(handle, format)
    else:
        # Get opener for compression if needed
        _open, mode = open_mode(path.suffix)
        # Return a generator of SeqRecords the way it was intended
        with _open(path, mode) as handle:
            yield from read_handle(handle, format)


def read_seq(handle: Union[str, StringIO], format: str = None) -> pd.DataFrame:
    """Reads a genomic file and converts it to a pandas DataFrame.

    Parameters
    ----------
    handle : str or StringIO
        Path to a genomic file or a StringIO object.
    format : str, optional
        File format. Default is fasta

    Returns
    -------
    pandas.DataFrame

    Examples
    --------
    >>> read_seq('file.fasta.gz')
    >>> read_seq('file.fasta.bz2')
    >>> read_seq('file.fa.gz')
    >>> read_seq('file.fq.gz')
    >>> read_seq('nosuffix_file', 'fasta')
    """
    seqrecords = read(handle, format)
    # Convert to pandas DataFrame - pulling out all the annotations and such as their own columns
    return BioDataFrame.from_seqrecords(seqrecords)
