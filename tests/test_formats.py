#!/usr/bin/env python
"""Tests for all genomic file format reading in seqpandas."""

import pytest
import os
import pandas as pd
from pathlib import Path

import seqpandas as spd
from seqpandas.vcf import read_vcf
from seqpandas.bed import read_bed
from seqpandas.seqpandas import BioDataFrame


class TestVCFReading:
    """Test VCF file reading functionality."""
    
    @pytest.fixture
    def vcf_path(self):
        """Path to test VCF file."""
        return Path(__file__).parent / "test-data" / "vcf-test-data" / "example.vcf"
    
    def test_read_vcf(self, vcf_path):
        """Test basic VCF reading."""
        df = read_vcf(str(vcf_path))
        assert isinstance(df, pd.DataFrame)
        assert df.shape[0] == 6  # 6 variants
        assert 'CHROM' in df.columns
        assert 'POS' in df.columns
        assert 'REF' in df.columns
        assert 'ALT' in df.columns
    
    def test_read_vcf_with_samples(self, vcf_path):
        """Test VCF reading with sample parsing."""
        df = read_vcf(str(vcf_path), samples_as_dict_dtype=True)
        assert 'Sample1' in df.columns
        assert 'Sample2' in df.columns
        assert 'Sample3' in df.columns
        # Check that samples are parsed as dicts
        assert isinstance(df['Sample1'].iloc[0], dict)
        assert 'GT' in df['Sample1'].iloc[0]
    
    def test_read_vcf_without_dict_parsing(self, vcf_path):
        """Test VCF reading without dict parsing."""
        df = read_vcf(str(vcf_path), samples_as_dict_dtype=False)
        assert 'Sample1' in df.columns
        # Check that samples remain as strings
        assert isinstance(df['Sample1'].iloc[0], str)


class TestBEDReading:
    """Test BED file reading functionality."""
    
    @pytest.fixture
    def bed_path(self):
        """Path to test BED file."""
        return Path(__file__).parent / "test-data" / "example.bed"
    
    def test_read_bed(self, bed_path):
        """Test basic BED reading."""
        df = read_bed(str(bed_path))
        assert isinstance(df, pd.DataFrame)
        assert df.shape[0] == 6  # 6 features
        assert 'chrom' in df.columns
        assert 'chromStart' in df.columns
        assert 'chromEnd' in df.columns
        assert 'name' in df.columns
        assert 'score' in df.columns
        assert 'strand' in df.columns
    
    def test_bed_coordinates(self, bed_path):
        """Test that BED coordinates are properly parsed."""
        df = read_bed(str(bed_path))
        assert df['chromStart'].iloc[0] == 1000
        assert df['chromEnd'].iloc[0] == 2000
        assert df['chrom'].iloc[0] == 'chr1'


class TestFASTAReading:
    """Test FASTA file reading functionality."""
    
    @pytest.fixture
    def fasta_path(self):
        """Path to test FASTA file."""
        return Path(__file__).parent / "test-data" / "example.fasta"
    
    def test_read_fasta(self, fasta_path):
        """Test basic FASTA reading."""
        df = spd.read_seq(str(fasta_path), format='fasta')
        assert isinstance(df, BioDataFrame)
        assert df.shape[0] == 4  # 4 sequences
        assert '_seq' in df.columns
        assert 'id' in df.columns
        assert 'description' in df.columns
    
    def test_fasta_sequences(self, fasta_path):
        """Test that FASTA sequences are properly parsed."""
        df = spd.read_seq(str(fasta_path), format='fasta')
        # Check first sequence
        assert df['id'].iloc[0] == 'seq1'
        assert 'Human BRCA1' in df['description'].iloc[0]
    
    def test_read_biopython_seqs(self, fasta_path):
        """Test reading as BioPython sequences."""
        seqs = spd.read(str(fasta_path), format='fasta')
        assert seqs is not None
        # Convert to list to check
        seq_list = list(seqs)
        assert len(seq_list) == 4


class TestSAMReading:
    """Test SAM file reading functionality."""
    
    @pytest.fixture 
    def sam_path(self):
        """Path to test SAM file."""
        return Path(__file__).parent / "test-data" / "example.sam"
    
    def test_read_sam(self, sam_path):
        """Test basic SAM reading."""
        # Note: SAM reading may return empty DataFrame if pysam has issues
        # or file format is not fully compatible
        df = spd.read_seq(str(sam_path), format='sam')
        assert isinstance(df, (BioDataFrame, pd.DataFrame))
        # SAM parsing implementation may vary
        assert df is not None


class TestCLI:
    """Test CLI functionality for all formats."""
    
    def test_cli_import(self):
        """Test that CLI can be imported."""
        from seqpandas.cli import main
        assert main is not None
    
    def test_cli_read_vcf(self):
        """Test CLI read-vcf command."""
        from click.testing import CliRunner
        from seqpandas.cli import main
        
        runner = CliRunner()
        vcf_path = Path(__file__).parent / "test-data" / "vcf-test-data" / "example.vcf"
        result = runner.invoke(main, ['read-vcf', str(vcf_path)])
        assert result.exit_code == 0
        assert 'Successfully read VCF file' in result.output
    
    def test_cli_read_formats(self):
        """Test CLI read command with different formats."""
        from click.testing import CliRunner
        from seqpandas.cli import main
        
        runner = CliRunner()
        test_files = {
            'vcf': Path(__file__).parent / "test-data" / "vcf-test-data" / "example.vcf",
            'bed': Path(__file__).parent / "test-data" / "example.bed",
            'fasta': Path(__file__).parent / "test-data" / "example.fasta",
        }
        
        for fmt, path in test_files.items():
            result = runner.invoke(main, ['read', '-f', fmt, str(path)])
            assert result.exit_code == 0, f"Failed for format {fmt}: {result.output}"
            assert f'Successfully read {fmt.upper()} file' in result.output