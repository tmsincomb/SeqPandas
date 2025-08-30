import multiprocessing

import numpy as np
from Bio import SeqIO
from pysam import AlignmentFile

# Import from seqpandas package structure
from .tools.pathing import pathing

# Try to import Cython extension, fall back to pure Python if not available
try:
    from .cython_numpy import to_array

    CYTHON_AVAILABLE = True
except (ImportError, ModuleNotFoundError):
    # Silently fall back to pure Python implementation
    # This is common in editable installs where the extension might not build
    CYTHON_AVAILABLE = False

    def to_array(inp):
        """Pure Python fallback for to_array function."""
        return np.array(inp, dtype=np.int16)

    # Only warn if explicitly requested
    import os

    if os.environ.get("SEQPANDAS_WARN_NO_CYTHON"):
        import warnings

        warnings.warn("Cython extension not available, using pure Python fallback for to_array")


def get_ref_seq(reference_file: str, format: str = "fasta"):
    reference_file = pathing(reference_file)
    seqrecords = list(SeqIO.parse(reference_file, format=format))
    if len(seqrecords) > 1:
        raise ValueError("ERROR :: Reference fasta must be a single complete sequence")
    return seqrecords[0].seq


def get_alignments(alignment_file, threads=multiprocessing.cpu_count()):
    return AlignmentFile(pathing(alignment_file), threads=threads)


class Pileup:
    """
    A class for creating and managing genomic pileup data from alignment files.

    This class processes BAM/SAM alignment files against a reference sequence
    to create a pileup matrix containing base counts, insertions, and deletions
    at each position.
    """

    basepile = {
        "A": [1, 0, 0, 0],
        "C": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "T": [0, 0, 0, 1],
    }
    baseindex = {
        "A": -4,
        "C": -3,
        "G": -2,
        "T": -1,
    }
    base_index = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
    }
    base_insert_index = {
        "A": 4,
        "C": 5,
        "G": 6,
        "T": 7,
    }

    def __repr__(self):
        return "Pileup 2D matrix - REF/deletion + ALT + ALT/insertion"

    def __init__(
        self,
        reference_file,
        alignment_file,
        minimum_alignment_coverage: float = 0.75,
        save_pileup_to_destination: str = None,
        use_saved_pileup: str = None,
        minimum_coverage: int = 30,
        heterozygous_threshold: float = 0.25,
        minimum_variant_radius: int = 15,
        overwrite_variant_calls: bool = False,
        **kwargs,
    ) -> None:
        # helper variables
        self.contig_start = None
        self.contig_end = None
        self.minimum_coverage = minimum_coverage
        self.heterozygous_threshold = heterozygous_threshold  # minimum difference a secondary alt needs to be for a variant to be called.
        self.minimum_variant_radius = minimum_variant_radius  # distance variants are allowed to be without being classified as complex
        self.mapping_coverages = []
        # alignment parameters
        self.minimum_alignment_coverage = minimum_alignment_coverage
        # core data
        self.reference_file = pathing(reference_file)
        self.alignment_file = pathing(alignment_file)
        self.ref_seq = get_ref_seq(reference_file)
        # to speed-up debugging
        if use_saved_pileup:
            try:
                print("=== Using Input Pileup ===")
                self.pileup = str(pathing(use_saved_pileup))
                # open it up if its a npz file
                if isinstance(self.pileup, str):
                    self.pileup = np.load(self.pileup)
                for alignment in get_alignments(self.alignment_file):
                    self.contig_start = alignment.reference_start
                    break
                if overwrite_variant_calls:
                    print("=== Overwriting Variant Calls === ")
                    self.overwrite_variant_calls()
            except (FileNotFoundError, ValueError, OSError) as e:
                print(f"=== Input Pileup Failed ({e}); Building Pileup From Scratch ===")
                self.pileup = self.get_contig_pileup()
        else:
            print("=== Building Pileup From Scratch ===")
            self.pileup = self.get_contig_pileup()
        if save_pileup_to_destination:
            self.save_pileup(save_pileup_to_destination)
        print("=== Pilup Complete ===")

    @staticmethod
    def _offset_sum(top, bot, left, right):
        top_len = len(top)
        bot_len = len(bot)
        if top.size == 0:
            return bot
        if right == 0 and left == 0:
            top += bot
        elif top_len <= left:
            padding = left - top_len
            return np.concatenate(
                (np.pad(top, (0, padding), "constant", constant_values=(0, 0)), bot)
            )
        elif right > 0:
            top[left:] += bot[: top_len - left]
            top = np.concatenate((top, bot[-right:]))
        else:
            top[left : bot_len + left] += bot[:]
        return top

    def _get_contig_pileup_fragment(self, alignment):
        align_pairs = alignment.get_aligned_pairs()[
            alignment.query_alignment_start : alignment.query_alignment_end
        ]
        contig_fragment = []
        for query_pos, ref_pos in align_pairs:
            if ref_pos is not None:
                last_ref_pos = ref_pos  # removes odd alignment edge-cases
            if ref_pos is not None and query_pos is not None:
                contig_fragment.extend(
                    self.basepile[alignment.query_sequence[query_pos]][:] + [0, 0, 0, 0]
                )
            elif (
                ref_pos is None and query_pos is not None and contig_fragment
            ):  # Insertion: query base with no reference position
                contig_fragment[self.baseindex[alignment.query_sequence[query_pos]]] += 1
            elif (
                ref_pos is not None and query_pos is None
            ):  # Deletion: reference position with no query base
                contig_fragment.extend([0, 0, 0, 0, 0, 0, 0, 0])
        return contig_fragment, last_ref_pos

    def _get_contig_pileup(self):
        pileup_len = 8
        alignmentfile = get_alignments(self.alignment_file)
        prev_reference_end = 0
        prev_contig_fragment = np.empty(1)
        for i, alignment in enumerate(alignmentfile):
            # If alignment less that 50% aligned, it's no good and will muddy up the pileup
            skip_base, total_aln_pos = 0, 0
            for code, adv in alignment.cigartuples:
                total_aln_pos += adv
                if code == 4:
                    skip_base += adv
            mapping_coverage = 1.0 - 1.0 * skip_base / (total_aln_pos + 1)
            if mapping_coverage < self.minimum_alignment_coverage:
                continue
            self.mapping_coverages.append(mapping_coverage)
            # initial start
            if i == 0 or not self.contig_start:
                self.contig_start = alignment.reference_start
                prev_reference_start = 0
                contig_pf, prev_reference_end = self._get_contig_pileup_fragment(alignment)
                prev_reference_end -= self.contig_start
                prev_contig_fragment = to_array(contig_pf)
                continue
            reference_start = alignment.reference_start - self.contig_start
            contig_pf, reference_end = self._get_contig_pileup_fragment(alignment)
            reference_end -= self.contig_start
            contig_fragment = to_array(contig_pf)
            # arrays need to to be the same size
            left = reference_start - prev_reference_start
            right = (
                reference_end - prev_reference_end
            )  # can be negative offset for additional padding on incoming alignment
            left *= pileup_len
            right *= pileup_len
            # padding arrays to be the same size
            merged_contig_fragment = self._offset_sum(
                prev_contig_fragment, contig_fragment, left, right
            )

            # alignments are sorted so the current reference slice is next
            prev_contig_fragment = merged_contig_fragment[
                (reference_start - prev_reference_start) * pileup_len :
            ]
            # previous reference start to current reference start slice for vstack
            fragment = merged_contig_fragment[
                : (reference_start - prev_reference_start) * pileup_len
            ]
            # update positions
            prev_reference_end = (
                reference_end if reference_end >= prev_reference_end else prev_reference_end
            )
            prev_reference_start = reference_start
            try:
                if fragment.size != 0:  # repeating alignments at same ref start position
                    yield fragment.reshape(-1, pileup_len)
            except ValueError as e:
                print(f"Error reshaping fragment at alignment {i}: {e}")
                print(
                    f"Fragment shape: {fragment.shape}, merged shape: {merged_contig_fragment.shape}"
                )
                raise RuntimeError(f"Failed to reshape pileup fragment at alignment {i}") from e
        self.contig_end = prev_reference_end + self.contig_start
        if prev_contig_fragment.size != 0:
            yield prev_contig_fragment.reshape(-1, pileup_len)
        alignmentfile.close()

    def _meets_filter(self, total_count: int, base_count: list, ref_base: str):
        """
        Determine if a position meets the criteria for variant calling.

        Args:
            total_count: Total number of reads covering this position
            base_count: List of counts for each base [A, C, G, T]
            ref_base: Reference base at this position

        Returns:
            int: 1 if position is a variant, 0 otherwise
        """
        base_count = zip("ACGT", base_count)
        base_count = sorted(base_count, key=lambda x: -x[1])
        if base_count[0][0] != ref_base:
            return 1
        p0 = 1.0 * base_count[0][1] / total_count
        p1 = 1.0 * base_count[1][1] / total_count
        if (
            p0 < 1.0 - self.heterozygous_threshold and p1 > self.heterozygous_threshold
        ):  # or base_count[0][0] != ref_base:
            return 1
        return 0

    #  TODO: There is a bottleneck somewhere in here. Memory issues are most likely the issue since time is linear with matrix depth & time allocating the actual memory was solved with cython.
    def get_contig_pileup(self):
        """
        Get pileup of reference+deletions, alt and alt+insertions for dynamic creation of tensors.

        Returns:
            csr_matrix: sparce column matrix to hold the complete segment of sequences in a tiny amount of memory.
        """
        contig_pileup = np.concatenate(list(self._get_contig_pileup()), axis=0)  # vstack: top+bot
        ref_sum = contig_pileup[:, :4].sum(axis=1)
        ref_seq_seg = self.ref_seq[self.contig_start :]
        ref_piles = []
        for pos, (nt, rs, base_count) in enumerate(zip(ref_seq_seg, ref_sum, contig_pileup[:, :4])):
            pile = [0] + self.basepile.get(nt, [0, 0, 0, 0])[:]
            try:  # TODO: give partial distributions for nts using http://www.bioinformatics.org/sms/iupac.html
                # if pos + self.contig_start + 1 == 14513747:
                #     print(pos, nt, rs, base_count)
                pile[self.base_index[nt] + 1] = rs
            except KeyError:
                # Handle unknown nucleotides (not A, C, G, T)
                ref_piles.append(pile)
                continue
            if rs >= self.minimum_coverage:
                pile[0] = self._meets_filter(rs, base_count, nt)
            ref_piles.append(pile)
        ref_pileup = np.array(ref_piles)
        return np.concatenate((ref_pileup, contig_pileup), axis=1)  # hstack: left+right

    def overwrite_variant_calls(self):
        ref_sum = self.pileup[:, 1:5].sum(axis=1)
        ref_seq_seg = self.ref_seq[self.contig_start :]
        ref_piles = []
        for pos, (nt, rs, base_count) in enumerate(zip(ref_seq_seg, ref_sum, self.pileup[:, 5:9])):
            variant_call = [0]
            if rs >= self.minimum_coverage:
                variant_call = [self._meets_filter(rs, base_count, nt)]
            ref_piles.append(variant_call)
        ref_pileup = np.array(ref_piles)
        self.pileup = np.concatenate((ref_pileup, self.pileup[:, 1:]), axis=1)  # hstack: left+right

    def save_pileup(self, destination):
        """
        Save scipy array to a npz file so we don't have to build the array from scratch.

        Args:
            destination ([type]): [description]
        """
        np.save(pathing(destination, new=True, overwrite=True), self.pileup)


if __name__ == "__main__":
    pileup = Pileup(
        reference_file="/home/tmsincomb/Desktop/SupplementaryDataset1/SIMULATED_SNP_CONTAINING_GENOMES/neisseria/FDAARGOS_205/FDAARGOS_205_simulated.fasta",
        alignment_file="/home/tmsincomb/Desktop/SupplementaryDataset1/SIMULATED_SNP_CONTAINING_GENOMES/neisseria/FDAARGOS_205/dwgsim-default.mapped.bam",
    )
    print(pileup.pileup[:10])
    # import IPython
    # IPython.embed()
    # print(pileup.pileup[:10])

    vcf = pathing("~/Dropbox/thesis/VariantNET/testing_data/chr21/chr21.vcf")
    fasta = pathing("~/Dropbox/thesis/VariantNET/testing_data/chr21/chr21.fa")
    bam = pathing(
        "~/Dropbox/thesis/VariantNET/testing_data/chr21/hg38.NA12878-WashU_chr21-14069662-46411975.bam"
    )
    bami = pathing(
        "~/Dropbox/thesis/VariantNET/testing_data/chr21/hg38.NA12878-WashU_chr21-14069662-46411975.bam.bai"
    )
    bed = pathing(
        "~/Dropbox/thesis/VariantNET/testing_data/chr21/CHROM21_v.3.3.2_highconf_noinconsistent.bed"
    )
    # p = Pileup(fasta, bam)
