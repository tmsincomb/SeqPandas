"""Cython optimized numpy operations for seqpandas."""

try:
    from .cython_np_array import to_array
    __all__ = ['to_array']
except ImportError:
    # Silently fall back to pure Python implementation
    # This is expected in editable installs or when Cython extension can't be built
    import numpy as np
    
    def to_array(inp):
        """Pure Python fallback for to_array function."""
        return np.array(inp, dtype=np.int16)
    
    __all__ = ['to_array']