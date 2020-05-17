import Bio
import pandas as pd
from pandas import DataFrame
import pandas
import numpy as np
from pandas.core.frame import *
import types
from typing import List


class SubclassedSeries(pd.Series):
    """ Pandas Series API to Inherit """
    @property
    def _constructor(self):
        return SubclassedSeries

    @property
    def _constructor_expanddim(self):
        return SubclassedDataFrame


class SubclassedDataFrame(pd.DataFrame):
    """ Pandas DataFrame to Inherit """
    @property
    def _constructor(self):
        return SubclassedDataFrame

    @property
    def _constructor_sliced(self):
        return SubclassedSeries


class BioDatabase(SubclassedDataFrame):
    """ Expanded Pandas DataFrame to handle BioPython SeqRecords generator or genomic file types """
    @classmethod
    def from_seqrecords(cls, seqrecords, index=None, exclude=None, columns=None,
                        coerce_float=False, nrows=None):
        """ Takes Biopython parsed output to convert to a proper DataFrame"""
        if isinstance(seqrecords, types.GeneratorType):
            data = cls.__normalize_seqrecords(seqrecords)
        else:
            data = seqrecords
        return cls.from_records(data, index=index, exclude=exclude, columns=columns,
                                coerce_float=coerce_float, nrows=nrows)

    def __normalize_seqrecords(seqrecords) -> List[dict]:
        """ Pull nested dictionaries into a single dictionary.

        Priority is given to the keys higher in the hierarchy.
        """
        records = []
        for seqrecord in SeqIO.to_dict(seqrecords).values():
            _records = []
            record = seqrecord.__dict__
            # If a more complicated format is used; features will be nested.
            features = record.pop('features') if record.get('features') else []
            for feature in features:
                _record = deepcopy(record)
                # Meta that make up the feature
                aspects = feature.__dict__
                # Qualifier dictionary inside each feature
                qualifiers = aspects.pop('qualifiers') if aspects.get('qualifiers') else {}
                # Add feature
                for aspect_key, aspect_value in aspects.items():
                    if aspect_key not in record:
                        _record[aspect_key] = aspect_value
                # Add qualifier
                for qualifier_key, qualifier_value in qualifiers.items():
                    _record = deepcopy(_record)
                    if qualifier_key not in _record:
                        _record[qualifier_key] = qualifier_value
                        _records += [_record]
                # If no qualifiers dump feature
                if not _records:
                    _records += [_record]
            # If no feature dump original seq record
            if not _records:
                _records += [record]
            records += _records

        return records


def read_seq(handle, format, alphabet=None):
    seqrecords = SeqIO.parse(handle, format=format, alphabet=alphabet)
    return BioDatabase.from_seqrecords(seqrecords)


pd.DataFrame = BioDatabase
pd.read_seq = read_seq
