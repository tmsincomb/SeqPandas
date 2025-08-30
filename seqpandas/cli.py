"""Console script for seqpandas."""

import sys
from pathlib import Path

import click

from seqpandas.vcf import read_vcf


@click.group()
def main():
    """SeqPandas CLI - Import genomic data into Pandas DataFrames."""
    pass


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--ignore-header", is_flag=True, help="Ignore VCF header and use default")
@click.option(
    "--samples-as-dict", is_flag=True, default=True, help="Parse sample columns as dictionaries"
)
@click.option("--output", "-o", type=click.Path(), help="Output CSV file (optional)")
def read_vcf_cmd(input_file, ignore_header, samples_as_dict, output):
    """Read a VCF file and display or save as CSV."""
    try:
        # Read the VCF file
        df = read_vcf(
            input_file, ignore_header=ignore_header, samples_as_dict_dtype=samples_as_dict
        )

        # Display basic info
        click.echo(f"Successfully read VCF file: {input_file}")
        click.echo(f"Shape: {df.shape[0]} variants, {df.shape[1]} columns")
        click.echo(f"Columns: {', '.join(df.columns)}")
        click.echo("\nFirst 5 variants:")
        click.echo(df.head().to_string())

        # Save to CSV if output specified
        if output:
            df.to_csv(output, index=False)
            click.echo(f"\nSaved to: {output}")

        return 0
    except Exception as e:
        click.echo(f"Error reading VCF file: {e}", err=True)
        return 1


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--format",
    "-f",
    type=click.Choice(["fasta", "sam", "vcf", "bed"]),
    required=True,
    help="File format",
)
def read(input_file, format):
    """Read various genomic file formats."""
    import seqpandas as spd

    try:
        if format == "vcf":
            df = spd.read_vcf(input_file)
        elif format == "bed":
            df = spd.read_bed(input_file)
        elif format in ["fasta", "sam"]:
            df = spd.read_seq(input_file, format=format)

        click.echo(f"Successfully read {format.upper()} file: {input_file}")
        click.echo(f"Shape: {df.shape}")
        click.echo("\nFirst 5 rows:")
        click.echo(df.head().to_string())
        return 0
    except Exception as e:
        click.echo(f"Error reading {format.upper()} file: {e}", err=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
