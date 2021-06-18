#!/usr/bin/env python3
from abc import ABC
from collections import OrderedDict
from copy import deepcopy  # copies all nested dictionaries
import gzip
from io import StringIO
import multiprocessing
from pathlib import Path
from types import GeneratorType
from typing import Union, Dict, List, Tuple
from sys import exit

import Bio
from Bio import SeqIO  # __init__ kicks in
import pandas
from pandas import DataFrame

try:
    from pysam import AlignmentFile, VariantFile
except:
    print('pysam disabled ::: conda install pysam for read_seq BAM functionality')

from .tools import pathing
# from pandas.core.frame import * # Needed if I want to dive even deeper.


class SubclassedSeries(pandas.Series):
    """ Pandas Series API to Inherit """

    @property
    def _constructor(self):
        return SubclassedSeries

    @property
    def _constructor_expanddim(self):
        return SubclassedDataFrame


class SubclassedDataFrame(pandas.DataFrame, ABC):
    """ Pandas DataFrame to Inherit """

    @property
    def _constructor(self):
        return SubclassedDataFrame

    @property
    def _constructor_sliced(self):
        return SubclassedSeries


# todo give a dataframe a method extension to become a tfrecord https://www.tensorflow.org/tutorials/load_data/tfrecord
class BioDataFrame(SubclassedDataFrame, ABC):
    """ Expanded Pandas DataFrame to handle BioPython SeqRecords generator or genomic file types """

    @classmethod
    def from_seqrecords(cls,
                        seqrecords: Union[GeneratorType, list],
                        index=None,
                        exclude=None,
                        columns=None,
                        coerce_float=False,
                        nrows=None) -> pandas.DataFrame:
        """ Takes Biopython parsed output to convert to a proper DataFrame

            See from_records for more details on the rest of the parameters
            # todo see if the rest should be included if we are going to customize it this way.

        :param seqrecords: Generator or list from BioPython universal output. All formats are the same output.
        :param index:
        :param exclude:
        :param columns:
        :param coerce_float:
        :param nrows:

        >>> from_seqrecords(Bio.SeqIO.parse('file.fasta', format='fasta'))
        """
        # Nomalize nested featues; will result in some redundant.
        # This won't be an issue since small seqrecord count usually means high
        # amount of features & vice versa .
        data = cls.__normalize_seqrecords(seqrecords)
        return cls.from_records(
            data,
            index=index,
            exclude=exclude,
            columns=columns,
            coerce_float=coerce_float,
            nrows=nrows
        )

    @staticmethod
    def __normalize_seqrecords(seqrecords: Union[GeneratorType, list]) -> List[dict]:
        """ Pull nested dictionaries into a single dictionary.

        Priority is given to the keys higher in the hierarchy.
        :param seqrecords: Generator from BioPython universal output. All formats are the same output.
        :returns: List of dictionaries with keys that were obtained along the way.

        >>> __normalize_seqrecords(Bio.SeqIO.parse('file.fasta', format='fasta'))
        """

        def update_record(_record, data, reference_keys):
            """ Quick update of dictionary without updating prexisting keys """
            for key, value in data.items():
                if key not in reference_keys:
                    _record[key] = value
            return _record

        seqrecords = SeqIO.to_dict(seqrecords).values()  # Only care for the actual records themselves
        records = []
        for seqrecord in seqrecords:
            _records = []
            # SeqIO Parse is a nested class
            record = seqrecord.__dict__
            # If a more complicated format is used; features will be nested.
            features = record.pop('features') if record.get('features') else []
            # Annotation dictionary inside each seqrecord
            annotations = record.pop('annotations') if record.get('annotations') else {}
            # Add each annotation as seperate column
            record = update_record(record, annotations, record.keys())
            for feature in features:
                _record = deepcopy(record)  # Needs to refresh for each feature
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

    # todo do I need this?
    def from_tf(self):
        pass


def read_bed(bedfile: str) -> pandas.DataFrame:
    """ Reads bed files without the pesky tbi index needed. 
        BED DOCS :: https://genome.ucsc.edu/FAQ/FAQformat.html
    """
    with open(bedfile, 'r') as f:
        header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
        data_start = 0
        for i, line in enumerate(f.read().split()):
            row = line.split('\t')
            if row[0] in ['browser', 'track']:
                data_start = i+1
                break
            if row[0] not in ['browser', 'track'] and i > 0:
                break
        return pandas.read_csv(bedfile, delimiter='\t', header=None, names=header, skiprows=data_start).fillna('')



def read_vcf(path: str, ignore_header: bool = False, samples_as_dict_dtype: bool = True,) -> pandas.DataFrame:
    """ VCF reader for Version 4.2 using https://samtools.github.io/hts-specs/VCFv4.2.pdf 

    Args:
        path (str): path to VCF file.
        ignore_header (bool, optional): Ignore header from VCF and use created one here. Defaults to False.
        samples_as_dict_dtype (bool, optional): convert samples str to dict using adjacent FORMAT column for the keys. 
            Defaults to True.

    Returns:
        pandas.DataFrame
    """
    # Required header; The next columns after info are FORMAT with each Sample{i} metadata trailing after
    header = OrderedDict([
        ('CHROM', str),  
        ('POS', int),  
        ('ID', str),                   
        ('REF', str),    
        ('ALT', str),  
        ('QUAL', str),                   
        ('FILTER', str), 
        ('INFO', str), 
        # FORMAT 
        # Sample1
        # Sample2
        # ...
        # SampleN
    ])
    input_header = None
    # Read VCF; pull header if possible
    with open(path, 'r') as f:
        lines = []
        for l in f:
            # hopefully remove EOF errors
            if not l: break
            # ignore meta header
            if l.startswith('##'): continue
            # use header if given
            if l.startswith('#'): 
                if not ignore_header:
                    input_header = OrderedDict([
                        (c, header.get(c, str)) # default to column being dtype str if its an extra
                        for c in l[1:-1].split('\t')  # start at 1 to remove the # in #CHROM; end at -1 to remove newline \n
                    ]) 
            else:
                lines.append(l)
    num_columns = len(lines[0].split('\t'))
    if input_header:
        num_header = len(input_header.keys())
        if num_header != num_columns:
            exit(f'ERROR :: You have a mismatch -> number of columns={num_columns} number of headers={num_header}')
        header = input_header
    else:            
        # If default header isn't enough, add them as Sample{i}.
        num_header = len(header.keys()) 
        extra_columns = num_columns - num_header
        if extra_columns == 1:
            exit(f'ERROR :: Column FORMAT was given with a training sample columns')
        elif extra_columns > 0:  # FORMAT column will not form unless it has sample columns with work with.
            header.update({'FORMAT': str}) # shared FORMAT of each trailing sample 
            header.update(OrderedDict([(f'Sample{i}', str) for i in range(1, extra_columns)]))
        else:
            pass  # No sample data to be used
    df = pandas.read_csv(
        StringIO(''.join(lines)),
        names=header.keys(),
        dtype=header,
        header=None, 
        sep='\t',
    )
    try:
        if not df.get('FORMAT'):
            df['FORMAT'] = 'GT'
            df['Sample1'] = ['0/1' if float(info['AF']) < 1 else '1/1' for info in df['INFO'].apply(lambda x: {v.split('=')[0]:v.split('=')[1] for v in x.split(';')})]
    except:
        pass  # attemp failed to create a useable genotype from allele frequency
    # try:
    if samples_as_dict_dtype:
        if not df['FORMAT'].empty:
            def uu(x):
                return dict(zip(x[0].split(':'), x[1].split(':')))
            for i in range(9, df.columns.shape[0]):
                df.iloc[:, i] = df.iloc[:, [8, i]].apply(lambda x: uu(x), axis=1)
    # except: 
    #     pass  # format column was improper
    return df


def read_seq(handle: Union[str, StringIO], format: str, threads=multiprocessing.cpu_count()) -> pandas.DataFrame:
    """ Read Bioinformatic file type

    :param str handle: str path of file to open.
    :param str format: Broad range of Bioinformatic formats ie fasta & genbank.
    :param int threads: Threads used for only the pysam parsers.

    >>> read_seq('file.fasta.gz', format='fasta')
    >>> read_seq('file.vcf', format='vcf')
    >>> read_seq('file.bcf', format='bcf')
    """
    format = format.lower().strip()
    # Only BioPython can handle StringIO for now.
    # if isinstance(handle, StringIO):
    #     seqrecords = SeqIO.read(handle, format=format)
    #     return BioDataFrame.from_seqrecords(seqrecords.__dict__)
    path = pathing(handle)
    # PySAM #
    try:
        if format == 'sam':
            samfile = AlignmentFile(path, 'r', threads=threads)
            return BioDataFrame([alignment.to_dict() for alignment in list(samfile)])
        if format == 'bam':
            samfile = AlignmentFile(path, 'rb', threads=threads)
            return BioDataFrame([alignment.to_dict() for alignment in list(samfile)])
        if format == 'cram':
            samfile = AlignmentFile(path, "rc", threads=threads)
            return BioDataFrame([alignment.to_dict() for alignment in list(samfile)])
        # TODO: VCF reader is bad anyway. archive this.
        # if format in ['vcf', 'bcf']:
        #     # WARNING! This shit will break if you don't have a header; use pandas.read_vcf if you want a simple df read.
        #     vcffile = VariantFile(path, threads=threads)
        #     header = vcffile.header.__str__().split('#')[-1].strip().split('\t')
        #     rows = [v.__str__().strip().split('\t') for v in vcffile]
        #     df = BioDataFrame(rows, columns=header)
        #     if not df['FORMAT'].empty:
        #         def uu(x):
        #             return dict(zip(x[0].split(':'), x[1].split(':')))
        #         for i in range(9, df.columns.shape[0]):
        #             df.iloc[:, i] = df.iloc[:, [8, i]].apply(lambda x: uu(x), axis=1)
        #     return df
    except:
        pass
    # BioPython #
    # If file is gzip compressed for BioPython
    if path.suffix == '.gz':
        with gzip.open(path, "rt") as handle:
            seqrecords = SeqIO.parse(handle, format=format)
            # need to use/return while I/O is open
            return BioDataFrame.from_seqrecords(seqrecords)
    # Uncompressed; will break if another compression is used besides gzip.
    seqrecords = SeqIO.parse(path, format=format)
    return BioDataFrame.from_seqrecords(seqrecords)


# ADDS to pandas module for seamless behavior #
pandas.DataFrame = BioDataFrame
pandas.read_seq = read_seq
pandas.read_vcf = read_vcf
pandas.read_bed = read_bed