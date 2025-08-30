project(
  'seqpandas',
  'c', 'cython',
  version: '0.0.2',
  license: 'MIT',
  meson_version: '>=1.2.0',
  default_options: [
    'buildtype=release',
    'b_ndebug=if-release',
    'c_std=c99',
  ],
)

# Python and NumPy dependencies
py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# NumPy dependency handling - robust approach for all scenarios
# 1. Try to find NumPy via dependency() first (most reliable)
# 2. Fall back to runtime detection
# 3. Allow build to continue even if NumPy not found (Cython module is optional)

# Method 1: Try to find NumPy via dependency system
np_dep = dependency('numpy', required: false)
numpy_found = false

if np_dep.found()
  # NumPy found via dependency system - most reliable
  inc_np = []  # Will be handled by np_dep
  numpy_found = true
  message('NumPy found via dependency system')
else
  # Method 2: Try to get NumPy include directory from Python
  incdir_numpy = run_command(py,
    ['-c', 'import numpy; print(numpy.get_include())'],
    check: false
  ).stdout().strip()
  
  if incdir_numpy != ''
    inc_np = include_directories(incdir_numpy)
    numpy_found = true
    message('NumPy include directory found via Python:', incdir_numpy)
  else
    # NumPy not available - Cython extension will be skipped or fail gracefully
    inc_np = []
    numpy_found = false
    warning('NumPy not found - Cython extension may not build')
  endif
endif

# Add subdirectories
subdir('seqpandas')