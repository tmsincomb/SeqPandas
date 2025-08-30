#!/usr/bin/env python
"""
Test tutorial functionality.
Validates that all tutorial code examples work correctly.
"""

import pytest
import pandas as pd
import numpy as np
import tempfile
import os

import seqpandas as spd
from seqpandas.alignment import (
    read_aligned_sam,
    show_alignment_block,
    alignment_to_frequency_matrix,
    generate_position_labels,
    get_alignment_frequencies,
    sequences_to_logo_format,
)


class TestTutorialFunctionality:
    """Test that tutorial examples work correctly."""
    
    @pytest.fixture
    def tutorial_alignment(self):
        """Example alignment from tutorial."""
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
    
    def test_read_aligned_sequences(self, sam_file, fasta_file):
        """Test Part 1: Reading aligned sequences."""
        # This should work even without reference file
        aligned_df = read_aligned_sam(
            alignment_file=sam_file,
            reference_file=None,  # Optional
            max_reads=10,
            show_insertions=True,
            gap_char="-",
            insertion_char="*",
        )
        
        assert isinstance(aligned_df, pd.DataFrame)
        assert len(aligned_df) <= 10
        assert 'aligned_sequence' in aligned_df.columns
        assert 'original_sequence' in aligned_df.columns
    
    def test_frequency_matrices(self, sam_file):
        """Test Part 1.5: Creating frequency matrices."""
        aligned_df = read_aligned_sam(sam_file, max_reads=10)
        
        if len(aligned_df) > 0:
            sequences = aligned_df["aligned_sequence"].tolist()
            
            # Create frequency matrix
            freq_matrix = alignment_to_frequency_matrix(sequences)
            assert isinstance(freq_matrix, pd.DataFrame)
            assert freq_matrix.shape[0] == len(sequences[0])
            
            # Generate position labels
            labels, indices = generate_position_labels(sequences[0])
            assert len(labels) == len(sequences[0])
            assert len(indices) == len(sequences[0])
            
            # Get alignment frequencies directly
            freq_df, position_labels = get_alignment_frequencies(
                sam_file, max_reads=10, normalize=True
            )
            assert isinstance(freq_df, pd.DataFrame)
            assert len(position_labels) == freq_df.shape[0]
            
            # Convert to logo format
            logo_df, logo_labels, variant_positions = sequences_to_logo_format(sequences)
            assert isinstance(logo_df, pd.DataFrame)
            assert len(logo_labels) == logo_df.shape[0]
            assert isinstance(variant_positions, list)
    
    @pytest.mark.skipif(
        not pytest.importorskip("logomaker", reason="logomaker not installed"),
        reason="logomaker required for logo tests"
    )
    def test_create_logos(self, tutorial_alignment):
        """Test Part 2: Creating sequence logos."""
        from seqpandas.logo import plot_logo_with_indels
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_logo:
            try:
                freq_df = plot_logo_with_indels(
                    alignment=tutorial_alignment,
                    output_file=tmp_logo.name,
                    title="Tutorial Test Logo",
                    color_scheme="NajafabadiEtAl2017",
                )
                
                assert isinstance(freq_df, pd.DataFrame)
                assert freq_df.shape[0] == len(tutorial_alignment[0])
                assert os.path.exists(tmp_logo.name)
                
                # Clean up CSV
                csv_file = tmp_logo.name.replace('.png', '.csv')
                if os.path.exists(csv_file):
                    os.unlink(csv_file)
            finally:
                if os.path.exists(tmp_logo.name):
                    os.unlink(tmp_logo.name)
    
    @pytest.mark.skipif(
        not pytest.importorskip("logomaker", reason="logomaker not installed"),
        reason="logomaker required for logo tests"
    )
    def test_sam_to_logo(self, sam_file):
        """Test Part 3: Creating logo from SAM alignments."""
        from seqpandas.logo import plot_logo_with_indels
        
        aligned_df = read_aligned_sam(sam_file, max_reads=10)
        
        if len(aligned_df) > 0 and "aligned_sequence" in aligned_df.columns:
            sequences = aligned_df["aligned_sequence"].tolist()
            sequences = [seq for seq in sequences if seq and len(seq) > 0][:10]
            
            if sequences:
                # Pad sequences to same length
                max_len = max(len(seq) for seq in sequences)
                sequences = [seq.ljust(max_len, "-") for seq in sequences]
                
                with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_logo:
                    try:
                        freq_df = plot_logo_with_indels(
                            alignment=sequences,
                            output_file=tmp_logo.name,
                            title="SAM Logo Test"
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
    
    def test_shannon_entropy_calculation(self, tutorial_alignment):
        """Test Part 4: Shannon entropy calculation."""
        freq_matrix = alignment_to_frequency_matrix(tutorial_alignment)
        
        def calculate_shannon_entropy(freq_matrix):
            entropies = []
            for idx, row in freq_matrix.iterrows():
                probs = row[row > 0]
                if len(probs) > 0:
                    entropy = -sum(probs * np.log2(probs))
                else:
                    entropy = 0
                entropies.append(entropy)
            return entropies
        
        entropies = calculate_shannon_entropy(freq_matrix)
        
        assert len(entropies) == freq_matrix.shape[0]
        assert all(0 <= e <= np.log2(5) for e in entropies)  # Max entropy for 5 bases
        
        mean_entropy = np.mean(entropies)
        max_entropy = np.max(entropies)
        
        assert 0 <= mean_entropy <= np.log2(5)
        assert 0 <= max_entropy <= np.log2(5)
    
    def test_consensus_sequence_generation(self, tutorial_alignment):
        """Test Part 4: Consensus sequence generation."""
        freq_matrix = alignment_to_frequency_matrix(tutorial_alignment)
        
        def generate_consensus_sequence(freq_matrix, threshold=0.5):
            consensus = []
            for idx, row in freq_matrix.iterrows():
                max_char = row.idxmax()
                max_freq = row.max()
                if max_freq >= threshold:
                    consensus.append(max_char if max_char != '-' else '')
                else:
                    consensus.append('N')
            return ''.join(consensus)
        
        consensus_seq = generate_consensus_sequence(freq_matrix, threshold=0.5)
        
        assert isinstance(consensus_seq, str)
        assert len(consensus_seq.replace('', '')) <= len(tutorial_alignment[0])
        
        # Most positions should have consensus (not 'N')
        n_count = consensus_seq.count('N')
        assert n_count < len(consensus_seq) / 2  # Less than half should be ambiguous
    
    def test_alignment_block_display(self, sam_file):
        """Test alignment block visualization."""
        aligned_df = read_aligned_sam(sam_file, max_reads=5)
        
        if len(aligned_df) > 0:
            alignment_block = show_alignment_block(aligned_df)
            
            assert isinstance(alignment_block, str)
            assert len(alignment_block) > 0
            
            # Check that read names are included
            for _, row in aligned_df.head(3).iterrows():
                assert row['read_name'] in alignment_block
    
    def test_complete_workflow(self, sam_file):
        """Test complete tutorial workflow."""
        # Step 1: Read aligned sequences
        aligned_df = read_aligned_sam(sam_file, max_reads=5)
        assert len(aligned_df) > 0
        
        # Step 2: Get sequences
        sequences = aligned_df['aligned_sequence'].tolist()
        assert all(isinstance(seq, str) for seq in sequences)
        
        # Step 3: Create frequency matrix
        freq_matrix = alignment_to_frequency_matrix(sequences)
        assert freq_matrix.shape[0] == len(sequences[0])
        
        # Step 4: Generate position labels
        labels, indices = generate_position_labels(sequences[0])
        assert len(labels) == len(sequences[0])
        
        # Step 5: Convert to logo format
        logo_df, logo_labels, variant_pos = sequences_to_logo_format(sequences)
        assert logo_df.shape[0] == len(logo_labels)
        
        # Step 6: Show alignment block
        block = show_alignment_block(aligned_df)
        assert isinstance(block, str)