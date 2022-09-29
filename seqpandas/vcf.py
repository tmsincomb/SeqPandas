from collections import OrderedDict
from datetime import datetime
from email import header
from io import StringIO

import pandas as pd

from seqpandas.seqpandas import BioDataFrame
from seqpandas import __version__


def read_vcf(
    path: str,
    ignore_header: bool = False,
    samples_as_dict_dtype: bool = True,
) -> pd.DataFrame:
    """VCF reader for Version 4.2 using https://samtools.github.io/hts-specs/VCFv4.2.pdf
    Args:
        path (str): path to VCF file.
        ignore_header (bool, optional): Ignore header from VCF and use created one here. Defaults to False.
        samples_as_dict_dtype (bool, optional): convert samples str to dict using adjacent FORMAT column for the keys.
            Defaults to True.
    Returns:
        pandas.DataFrame
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
                    input_header = OrderedDict(
                        [
                            (c, header.get(c, str))  # default to column being dtype str if its an extra
                            for c in line[1:-1].split(
                                "\t"
                            )  # start at 1 to remove the # in #CHROM; end at -1 to remove newline \n
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
        if not df.get("FORMAT"):
            df["FORMAT"] = "GT"
            df["Sample1"] = [
                "0/1" if float(info["AF"]) < 1 else "1/1"
                for info in df["INFO"].apply(lambda x: {v.split("=")[0]: v.split("=")[1] for v in x.split(";")})
            ]
    except KeyError:
        pass  # attempt failed to create a useable genotype from allele frequency; it's possible that it's valid.

    if samples_as_dict_dtype:
        if not df["FORMAT"].empty:

            def uu(x):
                return dict(zip(x[0].split(":"), x[1].split(":")))

            for i in range(9, df.columns.shape[0]):
                df.iloc[:, i] = df.iloc[:, [8, i]].apply(lambda x: uu(x), axis=1)

    return df


def to_vcf(df, reference:) -> str:
    header = f"""
    ##fileformat=VCFv4.3
    ##fileDate={datetime.today().isoformat().split('T')[0].replace('-', '')}
    ##source=SeqPandasV{__version__}
    """
    df = open_function[path.suffix](path, usecols=columns, dtype=str, keep_default_na=False)
    
    
class VCF(BioDataFrame):
    
    @property
    def E(cls, path, **kwargs):
        return cls(read_vcf(path, **kwargs))
    
    
    