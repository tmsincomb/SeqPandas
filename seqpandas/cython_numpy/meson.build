# Cython extension module for fast array operations
# This module provides optimized array conversion functions

# Check if we should build the Cython extension
# It's optional - the Python code has a fallback
build_cython = true

# Dependencies for the extension
deps = [py_dep]
incs = []

# Try to find NumPy - use multiple strategies
if np_dep.found()
  # NumPy found via dependency system from parent
  deps += np_dep
  message('NumPy found via dependency system for Cython module')
else
  # Try to get NumPy include directory from parent
  # Check if NumPy was found in parent meson.build
  if numpy_found and is_variable('inc_np')
    incs = inc_np
    message('Using NumPy include from parent meson.build')
  else
    # Final attempt: try to import NumPy directly
    incdir_numpy_local = run_command(py,
      ['-c', 'import numpy; print(numpy.get_include())'],
      check: false
    ).stdout().strip()
    
    if incdir_numpy_local != ''
      incs = include_directories(incdir_numpy_local)
      message('Found NumPy include via direct import:', incdir_numpy_local)
    else
      warning('NumPy headers not found - skipping Cython extension')
      build_cython = false
    endif
  endif
endif

# Build the Cython extension if possible
if build_cython
  py.extension_module(
    'cython_np_array',
    'cython_np_array.pyx',
    dependencies: deps,
    include_directories: incs,
    install: true,
    subdir: 'seqpandas/cython_numpy',
    c_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION']  # Suppress NumPy API warnings
  )
else
  message('Skipping Cython extension build - will use pure Python fallback')
endif

# Install the __init__.py file
py.install_sources(
  '__init__.py',
  subdir: 'seqpandas/cython_numpy'
)