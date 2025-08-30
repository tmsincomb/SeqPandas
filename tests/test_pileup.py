#!/usr/bin/env python
"""Tests for pileup module and Cython extensions."""

import pytest
import numpy as np
from pathlib import Path


class TestCythonExtension:
    """Test Cython numpy extension."""
    
    def test_to_array_import(self):
        """Test that to_array can be imported."""
        from seqpandas.cython_numpy import to_array
        assert callable(to_array)
    
    def test_to_array_functionality(self):
        """Test to_array converts list to numpy array correctly."""
        from seqpandas.cython_numpy import to_array
        
        test_list = [1, 2, 3, 4, 5]
        result = to_array(test_list)
        
        assert isinstance(result, np.ndarray)
        assert result.dtype == np.int16
        assert np.array_equal(result, np.array(test_list, dtype=np.int16))
    
    def test_to_array_empty_list(self):
        """Test to_array with empty list."""
        from seqpandas.cython_numpy import to_array
        
        result = to_array([])
        assert isinstance(result, np.ndarray)
        assert result.shape == (0,)
        assert result.dtype == np.int16
    
    def test_to_array_large_list(self):
        """Test to_array with large list."""
        from seqpandas.cython_numpy import to_array
        
        test_list = list(range(10000))
        result = to_array(test_list)
        
        assert result.shape == (10000,)
        assert result[0] == 0
        assert result[-1] == 9999


class TestPileupModule:
    """Test pileup module functionality."""
    
    def test_pileup_import(self):
        """Test that Pileup class can be imported."""
        try:
            from seqpandas import Pileup
            assert Pileup is not None
        except ImportError:
            # Direct import if scipy not available
            from seqpandas.pileup import Pileup
            assert Pileup is not None
    
    def test_pileup_utilities(self):
        """Test pileup utility functions."""
        from seqpandas.pileup import get_ref_seq, get_alignments
        
        assert callable(get_ref_seq)
        assert callable(get_alignments)
    
    def test_pileup_basepile(self):
        """Test Pileup class has correct base mappings."""
        from seqpandas.pileup import Pileup
        
        assert Pileup.basepile == {
            'A': [1, 0, 0, 0],
            'C': [0, 1, 0, 0],
            'G': [0, 0, 1, 0],
            'T': [0, 0, 0, 1],
        }
    
    def test_pileup_base_indices(self):
        """Test Pileup class has correct base indices."""
        from seqpandas.pileup import Pileup
        
        assert Pileup.base_index == {
            'A': 0,
            'C': 1,
            'G': 2,
            'T': 3,
        }
        
        assert Pileup.base_insert_index == {
            'A': 4,
            'C': 5,
            'G': 6,
            'T': 7,
        }