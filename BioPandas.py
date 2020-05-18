#!/usr/bin/env python3
from copy import deepcopy # copies all nested dictionaries
import gzip
from pathlib import Path
from types import GeneratorType
from typing import Union, Dict, List, Tuple

import Bio
from Bio import SeqIO # __init__ kicks in
from Bio.Alphabet import single_letter_alphabet # General sequence type
import pandas
from pandas import DataFrame
# from pandas.core.frame import * # Needed if I want to dive even deeper.


class SubclassedSeries(pandas.Series):
    """ Pandas Series API to Inherit """

    @property
    def _constructor(self):
        return SubclassedSeries

    @property
    def _constructor_expanddim(self):
        return SubclassedDataFrame


class SubclassedDataFrame(pandas.DataFrame):
    """ Pandas DataFrame to Inherit """

    @property
    def _constructor(self):
        return SubclassedDataFrame

    @property
    def _constructor_sliced(self):
        return SubclassedSeries


class DataFrame(SubclassedDataFrame):
    """ Expanded Pandas DataFrame to handle BioPython SeqRecords generator or genomic file types """

    @classmethod
    def from_seqrecords(cls, seqrecords: Union[GeneratorType, list], index=None, exclude=None,
                        columns=None, coerce_float=False, nrows=None) -> pandas.DataFrame:
        """ Takes Biopython parsed output to convert to a proper DataFrame

        :param seqrecords: Generator or list from BioPython universal output.
            All formats are the same output.

        >>> from_seqrecords(Bio.SeqIO.parse('file.fasta', format='fasta'))
        """
        # Nomalize nested featues; will result in some redundant.
        # This won't be an issue since small seqrecord count usually means high
        # amount of features & vice versa .
        data = cls.__normalize_seqrecords(seqrecords)
        return cls.from_records(data, index=index, exclude=exclude, columns=columns,
                                coerce_float=coerce_float, nrows=nrows)

    def __normalize_seqrecords(seqrecords: Union[GeneratorType, list]) -> List[dict]:
        """ Pull nested dictionaries into a single dictionary.

        Priority is given to the keys higher in the hierarchy.
        :param seqrecords: Generator from BioPython universal output.
            All formats are the same output.
        :returns: List of dictionaries with keys that were obtained along the way.

        >>> __normalize_seqrecords(Bio.SeqIO.parse('file.fasta', format='fasta'))
        """

        def update_record(_record, data, reference_keys):
            """ Quick update of dictionary without updating prexisting keys """
            for key, value in data.items():
                if key not in reference_keys:
                    _record[key] = value
            return _record

        records = []
        for seqrecord in SeqIO.to_dict(seqrecords).values():
            _records = []
            record = seqrecord.__dict__
            # If a more complicated format is used; features will be nested.
            features = record.pop('features') if record.get('features') else []
            # Annotation dictionary inside each seqrecord
            annotations = record.pop('annotations') if record.get('annotations') else {}
            # Add each annotation as seperate column
            record = update_record(record, annotations, record.keys())
            for feature in features:
                _record = deepcopy(record) # Needs to refresh for each feature
                # Meta that make up the feature
                aspects = feature.__dict__
                # Qualifier dictionary inside each feature
                qualifiers = aspects.pop('qualifiers') if aspects.get('qualifiers') else {}
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


def pathing(path: Union[str, Path], new: bool = False) -> Path:
    """ Guarantees correct expansion rules for pathing.

    :param Union[str, Path] path: path of folder or file you wish to expand.
    :param bool new: will check if distination exists if new
        (will check parent path regardless).
    :return: A pathlib.Path object.

    >>> pathing('~/Desktop/folderofgoodstuffs/')
    /home/user/Desktop/folderofgoodstuffs
    """
    path = Path(path)
    # Expand tilda shortened path or local path.
    if str(path)[0] == '~':
        path = path.expanduser()
    else:
        path = path.absolute() #
    # Making sure new paths don't exist while also making sure existing paths actually exist.
    if new:
        if not path.parent.exists():
            raise ValueError(f'ERROR ::: Parent directory of {path} does not exist.')
        if path.exists():
            raise ValueError(f'ERROR ::: {path} already exists!')
    else:
        if not path.exists():
            raise ValueError(f'ERROR ::: Path {path} does not exist.')
    return path


def read_seq(handle: str, format: str, alphabet: object = None) -> pandas.DataFrame:
    """ Read Bioinformatic file type

    :param str handle: str path of file to open.
    :param str format: Broad range of Bioinformatic formats ie fasta & genbank.
    :param object alphabet: Any number of Custom string object from BioPython with handy methods.
        [default: Bio.Alphabet.single_letter_alphabet]

    >>> read_seq('file.fasta.gz', format='fasta')
    """
    # Checks path validity
    path = pathing(handle)
    # If file is gzip compressed
    if path.suffix == '.gz':
        with gzip.open(handle, "rt") as handle:
            seqrecords = SeqIO.parse(handle, format=format, alphabet=alphabet)
            # need to use/return while I/O is open
            return DataFrame.from_seqrecords(seqrecords)
    # Uncompressed; will break if another compression is used.
    seqrecords = SeqIO.parse(handle, format=format, alphabet=alphabet)
    return DataFrame.from_seqrecords(seqrecords)


### Add to pandas module for seamless behavior ###
pandas.DataFrame = DataFrame
pandas.read_seq = read_seq
