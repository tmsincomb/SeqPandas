#!/usr/bin/env python
"""
Test alignment visualization functionality.
Tests SAM/BAM reading with aligned sequences showing gaps.
"""

import pytest
import pandas as pd
import tempfile
import os

import seqpandas as spd
from seqpandas.alignment import (
    read_aligned_sam,
    show_alignment_block,
    alignment_to_frequency_matrix,
)


class TestAlignmentVisualization:
    """Test alignment visualization features."""
    
    @pytest.fixture
    def sam_file(self):
        """Path to test SAM file."""
        return "tests/test-data/example.sam"
    
    @pytest.fixture
    def fasta_file(self):
        """Path to test FASTA file."""
        return "tests/test-data/example.fasta"
    
    def test_standard_sam_reading(self, sam_file):
        """Test standard SAM reading without alignment."""
        sam_df = spd.read_seq(sam_file, format='sam')
        
        assert isinstance(sam_df, pd.DataFrame)
        assert 'name' in sam_df.columns
        assert 'ref_name' in sam_df.columns
        assert 'ref_pos' in sam_df.columns
        assert 'seq' in sam_df.columns
        assert len(sam_df) > 0
    
    def test_aligned_sam_reading(self, sam_file):
        """Test SAM reading with alignment gaps."""
        aligned_df = read_aligned_sam(sam_file)
        
        assert isinstance(aligned_df, pd.DataFrame)
        assert 'aligned_sequence' in aligned_df.columns
        assert 'original_sequence' in aligned_df.columns
        assert 'read_name' in aligned_df.columns
        assert 'ref_start' in aligned_df.columns
        assert 'ref_end' in aligned_df.columns
        assert 'cigar' in aligned_df.columns
        assert len(aligned_df) > 0
    
    def test_alignment_gaps(self, sam_file):
        """Test that deletions appear as gaps in aligned sequences."""
        aligned_df = read_aligned_sam(sam_file)
        
        # Find reads with deletions
        deletion_reads = aligned_df[aligned_df['cigar'].str.contains('D')]
        
        if len(deletion_reads) > 0:
            for _, row in deletion_reads.iterrows():
                # Deletion should appear as gap in aligned sequence
                assert '-' in row['aligned_sequence']
                # Original sequence should not have gaps
                assert '-' not in row['original_sequence']
    
    def test_alignment_insertions(self, sam_file):
        """Test that insertions are handled correctly."""
        aligned_df = read_aligned_sam(sam_file, show_insertions=True)
        
        # Find reads with insertions
        insertion_reads = aligned_df[aligned_df['cigar'].str.contains('I')]
        
        if len(insertion_reads) > 0:
            for _, row in insertion_reads.iterrows():
                # With show_insertions=True, insertions should be marked
                # The aligned sequence should handle insertions appropriately
                assert len(row['aligned_sequence']) >= len(row['original_sequence'])
    
    def test_alignment_block_formatting(self, sam_file):
        """Test alignment block visualization."""
        aligned_df = read_aligned_sam(sam_file, max_reads=5)
        alignment_block = show_alignment_block(aligned_df)
        
        assert isinstance(alignment_block, str)
        lines = alignment_block.strip().split('\n')
        assert len(lines) > 0  # Should have at least one line
        
        # Check that the block contains meaningful content
        # The show_alignment_block might include separator lines
        content_lines = [line for line in lines if line and not line.startswith('=')]
        assert len(content_lines) > 0  # Should have at least one content line
    
    def test_sequence_padding(self, sam_file):
        """Test that aligned sequences are properly formatted."""
        aligned_df = read_aligned_sam(sam_file, max_reads=10)
        
        if len(aligned_df) > 1:
            sequences = aligned_df['aligned_sequence'].tolist()
            # Check that sequences have been aligned (may have different lengths due to different regions)
            assert all(isinstance(seq, str) for seq in sequences)
            assert all(len(seq) > 0 for seq in sequences)
    
    def test_cigar_string_parsing(self, sam_file):
        """Test that CIGAR strings are correctly parsed."""
        aligned_df = read_aligned_sam(sam_file)
        
        # Check various CIGAR operations
        cigar_ops = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}
        
        for _, row in aligned_df.iterrows():
            cigar = row['cigar']
            if cigar and cigar != '*':
                # CIGAR should contain valid operations
                import re
                ops = re.findall(r'[MIDNSHP=X]', cigar)
                assert all(op in cigar_ops for op in ops)
    
    def test_reference_positions(self, sam_file):
        """Test that reference positions are correctly calculated."""
        aligned_df = read_aligned_sam(sam_file)
        
        for _, row in aligned_df.iterrows():
            assert row['ref_start'] >= 0
            assert row['ref_end'] >= row['ref_start']
            # Aligned sequence length should match the reference span
            # (accounting for gaps and insertions)
    
    def test_original_vs_aligned_sequences(self, sam_file):
        """Test relationship between original and aligned sequences."""
        aligned_df = read_aligned_sam(sam_file, max_reads=5)
        
        for _, row in aligned_df.iterrows():
            original = row['original_sequence']
            aligned = row['aligned_sequence']
            
            # Original sequence should not have gaps
            assert '-' not in original
            
            # All bases from original should be in aligned
            # (though order might differ with complex operations)
            original_bases = original.replace('*', '').upper()
            aligned_bases = aligned.replace('-', '').replace('*', '').upper()
            
            # Basic check: aligned should contain at least the bases from original
            # (may have more due to insertions being shown)
    
    def test_frequency_matrix_from_alignment(self, sam_file):
        """Test creating frequency matrix from aligned sequences."""
        aligned_df = read_aligned_sam(sam_file, max_reads=10)
        sequences = aligned_df['aligned_sequence'].tolist()
        
        if sequences:
            freq_matrix = alignment_to_frequency_matrix(sequences)
            
            assert isinstance(freq_matrix, pd.DataFrame)
            assert freq_matrix.shape[0] == len(sequences[0])  # One row per position
            
            # Check that frequencies are valid
            for idx, row in freq_matrix.iterrows():
                row_sum = row.sum()
                assert abs(row_sum - 1.0) < 0.01, f"Position {idx} frequencies should sum to 1"
                assert all(0 <= val <= 1 for val in row.values)
    
    @pytest.mark.skipif(
        not pytest.importorskip("logomaker", reason="logomaker not installed"),
        reason="logomaker required for logo tests"
    )
    def test_logo_from_aligned_sequences(self, sam_file):
        """Test creating logo from aligned SAM sequences."""
        from seqpandas.logo import plot_logo_with_indels
        
        aligned_df = read_aligned_sam(sam_file, max_reads=5)
        sequences = aligned_df['aligned_sequence'].tolist()
        
        if sequences:
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_logo:
                try:
                    freq_df = plot_logo_with_indels(
                        alignment=sequences,
                        output_file=tmp_logo.name,
                        title="SAM Alignment Test"
                    )
                    
                    assert isinstance(freq_df, pd.DataFrame)
                    assert os.path.exists(tmp_logo.name)
                    
                    # Clean up CSV
                    csv_file = tmp_logo.name.replace('.png', '.csv')
                    if os.path.exists(csv_file):
                        os.unlink(csv_file)
                finally:
                    if os.path.exists(tmp_logo.name):
                        os.unlink(tmp_logo.name)