import pandas as pd


def read_bed(bedfile: str) -> pd.DataFrame:
    """Reads bed files without the pesky tbi index needed.
    BED DOCS :: https://genome.ucsc.edu/FAQ/FAQformat.html
    """
    with open(bedfile, "r") as f:
        header = [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "blockCount",
            "blockSizes",
            "blockStarts",
        ]
        data_start = 0
        for i, line in enumerate(f.read().split()):
            row = line.split("\t")
            if row[0] in ["browser", "track"]:
                data_start = i + 1
                break
            if row[0] not in ["browser", "track"] and i > 0:
                break
        return pd.read_csv(bedfile, delimiter="\t", header=None, names=header, skiprows=data_start).fillna("")
