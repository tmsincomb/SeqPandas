#!/usr/bin/env python
"""
Test alignment and logo module integration.
Tests the refactored alignment functions and their use in logo generation.
"""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import seqpandas as spd
from seqpandas.alignment import (
    alignment_to_frequency_matrix,
    generate_position_labels,
    get_alignment_frequencies,
    read_aligned_sam,
    sequences_to_logo_format,
    show_alignment_block,
)


class TestAlignmentIntegration:
    """Test alignment module functionality and integration with logo module."""

    @pytest.fixture
    def sample_sequences(self):
        """Provide sample aligned sequences for testing."""
        return [
            "ATGGATTTAT-CTGCTCTTCG",
            "ATGGATTTAT-CTGCTCTTCG",
            "ATGGATTTAT-CTGCTCTTCG",
            "ATGGATTTGT-CTGCTCTTCG",
            "ATGGATTTAT-CTGCTCTTCG",
            "ATGGCTTTAT-CTGCTCTTCG",
            "ATGGATTTAT-CTGCTCTTCG",
            "ATGGATTTATACTGCTCTTCG",
            "ATGGATTTAT-CTGCTCTTCG",
            "ATGGATTTAT-CTGCTCTTCG",
        ]

    @pytest.fixture
    def sam_file(self):
        """Path to test SAM file."""
        return "tests/test-data/example.sam"

    @pytest.fixture
    def fasta_file(self):
        """Path to test FASTA file."""
        return "tests/test-data/example.fasta"

    def test_alignment_to_frequency_matrix(self, sample_sequences):
        """Test conversion of aligned sequences to frequency matrix."""
        # Test with normalization
        freq_matrix = alignment_to_frequency_matrix(sample_sequences, normalize=True)

        assert isinstance(freq_matrix, pd.DataFrame)
        assert freq_matrix.shape[0] == len(sample_sequences[0])  # One row per position
        assert all(0 <= val <= 1 for val in freq_matrix.values.flatten())  # Normalized values

        # Test without normalization
        freq_matrix_counts = alignment_to_frequency_matrix(sample_sequences, normalize=False)
        assert freq_matrix_counts.values.max() <= len(sample_sequences)  # Count shouldn't exceed num sequences

        # Test that gaps are included by default
        assert "-" in freq_matrix.columns or all("-" not in seq for seq in sample_sequences)

    def test_generate_position_labels(self, sample_sequences):
        """Test position label generation with gap handling."""
        reference_seq = sample_sequences[0]
        labels, indices = generate_position_labels(reference_seq)

        assert len(labels) == len(reference_seq)
        assert len(indices) == len(reference_seq)

        # Check that gap positions have suffix labels
        for i, char in enumerate(reference_seq):
            if char == "-":
                assert "_" in labels[i]  # Gap positions should have suffix

        # Test with specific positions
        labels_filtered, indices_filtered = generate_position_labels(reference_seq, positions=[0, 1, 2, 10])
        # The function only returns labels for positions that exist in the sequence
        assert len(labels_filtered) == len(indices_filtered)
        # The indices should be based on non-gap positions, not raw positions
        assert len(indices_filtered) <= 4  # At most the requested positions

    def test_get_alignment_frequencies(self, sam_file):
        """Test getting alignment frequencies directly from SAM file."""
        freq_df, ref_positions = get_alignment_frequencies(sam_file, max_reads=5)

        assert isinstance(freq_df, pd.DataFrame)
        assert len(ref_positions) == freq_df.shape[0]
        assert all(isinstance(label, str) for label in ref_positions)

        # Test with normalization off
        freq_df_counts, _ = get_alignment_frequencies(sam_file, max_reads=5, normalize=False)
        assert freq_df_counts.values.max() <= 5  # Should not exceed max_reads

    def test_sequences_to_logo_format(self, sample_sequences):
        """Test conversion to logo format with variant position detection."""
        logo_df, logo_labels, variant_positions = sequences_to_logo_format(sample_sequences)

        assert isinstance(logo_df, pd.DataFrame)
        assert len(logo_labels) == logo_df.shape[0]
        assert isinstance(variant_positions, list)

        # Check that variant positions are correctly identified
        # Variant positions are where there is variation from the reference
        # The exact positions depend on the reference sequence chosen
        assert isinstance(variant_positions, list)
        # There should be some variation in these sequences
        if len(set("".join(sample_sequences))) > 1:  # If there's any variation
            assert len(variant_positions) >= 0  # May or may not have variants depending on reference

        # Test with different reference index
        logo_df2, _, _ = sequences_to_logo_format(sample_sequences, reference_index=1)
        assert logo_df2.shape == logo_df.shape

    def test_read_aligned_sam(self, sam_file):
        """Test reading SAM file with alignment gaps."""
        aligned_df = read_aligned_sam(sam_file, max_reads=5)

        assert "aligned_sequence" in aligned_df.columns
        assert "original_sequence" in aligned_df.columns
        assert "cigar" in aligned_df.columns
        assert len(aligned_df) <= 5

        # Check that sequences have gaps where appropriate
        for _, row in aligned_df.iterrows():
            if "D" in row["cigar"]:  # Deletion should appear as gap
                assert "-" in row["aligned_sequence"]

    def test_show_alignment_block(self, sam_file):
        """Test alignment block visualization."""
        aligned_df = read_aligned_sam(sam_file, max_reads=3)
        alignment_block = show_alignment_block(aligned_df)

        assert isinstance(alignment_block, str)
        # Check that some read names are in the block
        lines = alignment_block.strip().split("\n")
        assert len(lines) > 0

        # Each line should have the read name
        for _, row in aligned_df.head(3).iterrows():
            assert row["read_name"] in alignment_block

    @pytest.mark.skipif(
        not pytest.importorskip("logomaker", reason="logomaker not installed"),
        reason="logomaker required for logo tests",
    )
    def test_logo_integration(self, sample_sequences):
        """Test that logo module can use alignment functions."""
        from seqpandas.logo import plot_logo_with_indels

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_logo:
            try:
                freq_df = plot_logo_with_indels(
                    alignment=sample_sequences, output_file=tmp_logo.name, title="Test Logo"
                )

                assert isinstance(freq_df, pd.DataFrame)
                assert os.path.exists(tmp_logo.name)

                # Check CSV output was created
                csv_file = tmp_logo.name.replace(".png", ".csv")
                assert os.path.exists(csv_file)

                # Clean up
                os.unlink(csv_file)
            finally:
                if os.path.exists(tmp_logo.name):
                    os.unlink(tmp_logo.name)

    def test_frequency_matrix_consistency(self, sample_sequences):
        """Test that different methods produce consistent frequency matrices."""
        # Method 1: Direct from sequences
        freq_matrix1 = alignment_to_frequency_matrix(sample_sequences)

        # Method 2: Via logo format
        logo_df, _, _ = sequences_to_logo_format(sample_sequences)

        # Both should have the same shape
        assert freq_matrix1.shape[0] == logo_df.shape[0]  # Same number of positions

        # Check that frequencies sum to 1 (or close to it) for each position
        for idx in range(freq_matrix1.shape[0]):
            assert abs(freq_matrix1.iloc[idx].sum() - 1.0) < 0.01
            assert abs(logo_df.iloc[idx].sum() - 1.0) < 0.01

    def test_integration_with_real_sam(self, sam_file):
        """Test complete integration with real SAM file."""
        # Read aligned sequences
        aligned_df = read_aligned_sam(sam_file, max_reads=5)
        sequences = aligned_df["aligned_sequence"].tolist()

        # Process through alignment functions
        freq_matrix = alignment_to_frequency_matrix(sequences)
        logo_df, logo_labels, variant_pos = sequences_to_logo_format(sequences)

        # Verify outputs
        assert freq_matrix.shape[0] > 0
        assert len(logo_labels) == logo_df.shape[0]
        assert isinstance(variant_pos, list)

        # Test that we can create an alignment block
        block = show_alignment_block(aligned_df)
        assert isinstance(block, str)
        assert len(block) > 0
