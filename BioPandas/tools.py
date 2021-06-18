"""
Home for string tools and any other misc tool until its large enough for it's own file
"""
from subprocess import check_output, CalledProcessError
from typing import Union, Tuple
from multiprocessing import cpu_count
from pathlib import Path
import os
import sys


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def pathing(path: str, new: bool = False, overwrite: bool = True) -> Path:
    """ Guarantees correct expansion rules for pathing.

    :param Union[str, Path] path: path of folder or file you wish to expand.
    :param bool new: will check if distination exists if new  (will check parent path regardless).
    :return: A pathlib.Path object.

    >>> pathing('~/Desktop/folderofgoodstuffs/')
    /home/user/Desktop/folderofgoodstuffs
    """
    path = Path(path)
    # Expand shortened path
    if str(path)[0] == '~':
        path = path.expanduser()
    # Exand local path
    if str(path)[0] == '.':
        path = path.resolve()
    else:
        path = path.absolute()
    # Making sure new paths don't exist while also making sure existing paths actually exist.
    if new:
        if not path.parent.exists():
            raise ValueError(f'ERROR ::: Parent directory of {path} does not exist.')
        if path.exists() and not overwrite:
            raise ValueError(f'ERROR ::: {path} already exists!')
    else:
        if not path.exists():
            raise ValueError(f'ERROR ::: Path {path} does not exist.')
    return path


def mockreads(fasta: str, verbose=True) -> Tuple[str, str]:
    fasta = pathing(fasta)
    forward_read_path, backward_read_path = fasta.with_suffix('.read1.fq'), fasta.with_suffix('.read2.fq')
    cmd = f"wgsim -e 0 -r 0 -R 0 -X 0 -A 0 {fasta} {forward_read_path} {backward_read_path}"
    if verbose:
        print(cmd)
    _ = check_output(cmd, shell=True)
    return forward_read_path, backward_read_path


def bcftools_vcf(ref_file: str, bam_file: str, output_folder: str, output_prefix: str) -> str:
    ref_file = pathing(ref_file)
    bam_file = pathing(bam_file)
    output_folder = pathing(output_folder)
    output = output_folder / (output_prefix + '.bcftools.vcf')
    cmd = (
        f"bcftools mpileup -a AD -Ou --threads {cpu_count()} -f {ref_file} {bam_file} "
        f"| bcftools call -Ou -mv "
        f"| bcftools filter -s LowQual -e '%QUAL<20' > {output}"
    )
    try:
        _ = check_output(cmd, shell=True)
    except Exception as e:
        print(e)
    return output


def dwgsim(ref_file: str, output_folder: str, output_prefix: str, verbose: bool = True, **options) -> None:
    """
    dwgsim wrapper for simulating hom, het withing reads given a fasta file.

    Args:
        ref_file (str): reference fasta file.
        output_folder (str): all dwgsim outputs go into this folder
        output_prefix (str): output prefix for dwgsim to create simulated reads and mutation files.
    """
    ref_file = pathing(ref_file)
    output_folder = pathing(output_folder)
    hard_coded_options = {
        # '-r': .005,
        # '-X': 0,
        # '-R': 0,
        # '-1': 150,
        # '-2': 150,
        # '-C': 50,
        # '-H': True 
    }
    options = ' '.join([f'{k} {v}' if not isinstance(v, bool) else f'{k}' for k, v in {**hard_coded_options, **options}.items()])
    # unix "cd" simulator
    with cd(output_folder):
        cmd = f'dwgsim {options} {ref_file} {output_prefix}'
        if verbose: print(cmd)
        try:
            output = check_output(cmd, shell=True)
        except CalledProcessError as e:
            sys.exit(e.output)
