from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime
from email import header
from io import StringIO
from pathlib import Path

import pandas as pd

from .seqpandas import BioDataFrame

__version__ = "0.0.2"  # Import directly to avoid circular import

@dataclass
class IO:
    """ VCF I/O for VCF Version 4.3 found at src https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf
    
    META
        Specifics on how the VCF was generated.
    HEADER
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 Sample3
    BODY
        Variant call row following format 
    """
    meta: str = None
    
    # Required header; The next columns after info are FORMAT with each Sample{i} metadata trailing after
    default_header = OrderedDict(
        [
            ("CHROM", str),
            ("POS", int),
            ("ID", str),
            ("REF", str),
            ("ALT", str),
            ("QUAL", str),
            ("FILTER", str),
            ("INFO", str),
            # FORMAT
            # Sample1
            # Sample2
            # ...
            # SampleN
        ]
    )
    
    def get_vcf_meta(self, path: Path | str) -> list[str]:
        with open(path, "r") as infile:
            lines = []
            for line in infile:
                # hopefully removes EOF errors
                if not line:
                    break
                # ignore meta header
                if line.startswith("##"):
                    continue
                # use header if given
                if line.startswith("#"):
                    if ignore_header is False:
                        # start at 1 to remove the # in #CHROM; end at -1 to remove newline \n
                        column_names = line[1:-1].split("\t") 
                        column_name_type_mappings = OrderedDict(
                            [
                                (column_name, header.get(column_name, str))  # default to column being dtype str if its an extra
                                for column_name in column_names
                            ]
                        )
                else:
                    lines.append(line)
            num_columns = len(lines[0].split("\t"))
        
        if input_header:
            num_header = len(input_header.keys())
            if num_header != num_columns:
                exit(f"ERROR :: You have a mismatch -> number of columns={num_columns} number of headers={num_header}")
            header = input_header
        else:
            # If default header isn't enough, add them as Sample{i}.
            num_header = len(header.keys())
            extra_columns = num_columns - num_header
            if extra_columns == 1:
                exit("ERROR :: Column FORMAT was given with a training sample columns")
            elif extra_columns > 0:  # FORMAT column will not form unless it has sample columns with work with.
                header.update({"FORMAT": str})  # shared FORMAT of each trailing sample
                header.update(OrderedDict([(f"Sample{i}", str) for i in range(1, extra_columns)]))
            else:
                pass  # No sample data to be used

def read_vcf(
    path: str,
    ignore_header: bool = False,
    samples_as_dict_dtype: bool = True,
) -> pd.DataFrame:
    """VCF reader for Version 4.2 using https://samtools.github.io/hts-specs/VCFv4.2.pdf
    
    Parameters
    ----------
    path : str
        Path to a VCF file.
    ignore_header : bool, default False
        Ignore header from VCF and use created one here. Defaults to False.
        VCF are hard enough to parse; headers are rarely done right so we need to allow a way to ignore them.
        
    
    Returns
    -------
        pd.DataFrame
    """
    # Required header; The next columns after info are FORMAT with each Sample{i} metadata trailing after
    header = OrderedDict(
        [
            ("CHROM", str),
            ("POS", int),
            ("ID", str),
            ("REF", str),
            ("ALT", str),
            ("QUAL", str),
            ("FILTER", str),
            ("INFO", str),
            # FORMAT
            # Sample1
            # Sample2
            # ...
            # SampleN
        ]
    )
    input_header = None
    # Read VCF; pull header if possible
    with open(path, "r") as infile:
        lines = []
        for line in infile:
            # hopefully remove EOF errors
            if not line:
                break
            # ignore meta header
            if line.startswith("##"):
                continue
            # use header if given
            if line.startswith("#"):
                if not ignore_header:
                    # start at 1 to remove the # in #CHROM; end at -1 to remove newline \n
                    column_names = line[1:-1].split("\t") 
                    input_header = OrderedDict(
                        [
                            (column_name, header.get(column_name, str))  # default to column being dtype str if its an extra
                            for column_name in column_names
                        ]
                    )
            else:
                lines.append(line)
    num_columns = len(lines[0].split("\t"))
    if input_header:
        num_header = len(input_header.keys())
        if num_header != num_columns:
            exit(f"ERROR :: You have a mismatch -> number of columns={num_columns} number of headers={num_header}")
        header = input_header
    else:
        # If default header isn't enough, add them as Sample{i}.
        num_header = len(header.keys())
        extra_columns = num_columns - num_header
        if extra_columns == 1:
            exit("ERROR :: Column FORMAT was given with a training sample columns")
        elif extra_columns > 0:  # FORMAT column will not form unless it has sample columns with work with.
            header.update({"FORMAT": str})  # shared FORMAT of each trailing sample
            header.update(OrderedDict([(f"Sample{i}", str) for i in range(1, extra_columns)]))
        else:
            pass  # No sample data to be used
    df = pd.read_csv(
        StringIO("".join(lines)),
        names=header.keys(),
        dtype=header,
        header=None,
        sep="\t",
    )
    try:
        if "FORMAT" not in df.columns:
            df["FORMAT"] = "GT"
            df["Sample1"] = [
                "0/1" if float(info.get("AF", "1")) < 1 else "1/1"
                for info in df["INFO"].apply(lambda x: {v.split("=")[0]: v.split("=")[1] for v in x.split(";") if "=" in v})
            ]
    except (KeyError, ValueError):
        pass  # attempt failed to create a useable genotype from allele frequency; it's possible that it's valid.

    if samples_as_dict_dtype and "FORMAT" in df.columns:
        if not df["FORMAT"].empty:

            def uu(x):
                return dict(zip(x.iloc[0].split(":"), x.iloc[1].split(":")))

            for i in range(9, df.columns.shape[0]):
                df.iloc[:, i] = df.iloc[:, [8, i]].apply(lambda x: uu(x), axis=1)

    return df


def to_vcf(df, reference: str) -> str:
    from datetime import datetime
    header = f"""##fileformat=VCFv4.3
##fileDate={datetime.today().isoformat().split('T')[0].replace('-', '')}
##source=SeqPandasV{__version__}
##reference={reference}
"""
    # TODO: Implement VCF writing logic
    return header
    
    
class VCF(BioDataFrame):
    
    @classmethod
    def from_vcf(cls, path, **kwargs):
        return cls(read_vcf(path, **kwargs))
    
    
    